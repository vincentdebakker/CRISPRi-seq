import numpy as np
from numba import njit
import time
import multiprocessing
import datetime
from numba import types
from numba.typed import Dict
import os
import sys

# Type-expressions are currently not supported inside jit functions. 
# Leaving as global
int_array = types.int32[:]
inner_dict_type = types.DictType(types.unicode_type, types.int32[:])

@njit
def numba_dict(tipo):
    
    """ Creates a numba type dictionary """
    
    dictionary = Dict.empty(
        key_type=types.unicode_type,
        value_type=tipo)
    return dictionary

@njit
def sgrna_all_vs_all_comparator_inner(guide,all_combos,parameters):
    
    """ Runs the inner loop of the all vs all sgRNA comparison.
    Compares one sgRNA to all the PAM sites. 
    Returns the mismatch matrix"""
    
    result_dict = numba_dict(int_array)
    for b in all_combos: 
        result = binary_subtract(guide,all_combos[b],parameters) 
        if result is not None:
            result_dict[b]=result

    return result_dict

@njit
def sgrna_all_vs_all_comparator_outer(chunk,all_combos,parameters,i):
    
    """ Runs the outer loop of the all vs all sgRNA comparison.
    Sends individually the sgRNAs for the inner loop.
    Returns the final mismatch matrix"""
    
    compiled = numba_dict(inner_dict_type)
    for a in chunk:
        compiled[a] = sgrna_all_vs_all_comparator_inner(chunk[a],all_combos,parameters)
    return compiled

@njit
def binary_subtract(array1,array2,parameters):
    
    """ subtracts the binary matrix corresponding to one sgRNA, 
    to the comparing one. All the user inputed parameters are run, one by one,
    in an additive manner. if more mismatches than the allowed are detected,
    the function returns None, and the iteration for this combination is stopped"""

    for length, mismatch in parameters:
        sub = np.ones((len(array1[:length])), dtype=np.int32)
        miss=0
        for i,(arr1,arr2) in enumerate(zip(array1[:length],array2[:length])):
            if arr1-arr2 == 0:
                sub[i] = 0
            else:
                miss += 1           
            if miss>mismatch:
                return #if any of the rules are broken, ignore the entire thing

    return(sub)

def file_parser(file,nb_dict):
    
    """ Parses the input .csv files into dictionaries. Converts all DNA
    sequences to their respective binary array forms """
    
    def str_to_binary(text):
        text = list(map(bin,bytearray(text,'utf8')))
        for i, letter in enumerate(text):
            text[i] = int(letter[2:])
        return text
    
    if nb_dict:
        container = numba_dict(int_array)
    else:
        container = {}
    
    with open(file) as current: 
        for line in current:
            line = line[:-1].split(",")
            name = line[0]
            sequence = line[1].upper()
            sequence = sequence.replace(" ", "")
            byte_list = str_to_binary(sequence)
            if name not in container:
                container[name] = np.array((byte_list))
                
    return container

def cpu_count():
    
    """Counts the the ammount of cpu and 
    initializes the pool for multiprocessing """
    
    cpu = multiprocessing.cpu_count()
    if cpu >= 2:
        cpu -= 1
    pool = multiprocessing.Pool(processes = cpu)
    return pool,cpu

def processing_decision(sgrna,parameters,directory,names):
    
    """Runs a small scale simulation to infer total running time.
    Prepares the files for either multiprocessing or single processing"""
    
    pool,cpu = cpu_count()
    all_combos = file_parser(directory+names[1],True)
     
    #extrapolate total tunning time from the first 30 entries
    if len(sgrna)>=30:
        time_estimate(sgrna,all_combos,parameters,cpu)
    
    #single processed
    if len(sgrna) <= cpu:
        result=sgrna_all_vs_all_comparator_outer(sgrna,all_combos,parameters,1)
        write(result,directory+names[2].format(1))
    
    #multi processed
    else:
        chunked = dict_split(sgrna, cpu)
        out_paths=multi(chunked,parameters,pool,directory,names)
        out = directory+names[2].format("end")
        write_compiled(out_paths,out)

def time_estimate(sgrna,all_combos,parameters,cpu):
    
    """estimates the minimal total running time from the running time of a 
    small subset (30) of sgRNAs """
    
    #create estimate dictionary
    estimate = numba_dict(int_array)
    
    for key in sgrna:
        if len(estimate) < 30:
            estimate[key] = sgrna[key]
        else:
            break
        
    start = time.time()
    sgrna_all_vs_all_comparator_outer(estimate,all_combos,parameters,"test")
    end = time.time()-start
    
    estimate = end * (len(sgrna))/30/cpu
    end_date = datetime.datetime.now()+datetime.timedelta(seconds=estimate)
    
    print(f"\nThe estimated end time should be no sooner than {end_date.strftime('%c')}")

def multi_numba_compiler(chunk,parameters,i,directory,names):
    
    """Initializes the numba dictionaries with all the sgRNAs, and PAM sites.
    Because numba objects cannot be pickled at the moment, 
    all numba processing needs to be done after child process is created"""
    
    all_combos = file_parser(directory+names[1],True)
    chunk_numba = numba_dict(int_array)
    
    for chunky in chunk:
        chunk_numba[chunky] = chunk[chunky]
    
    result = sgrna_all_vs_all_comparator_outer(chunk_numba,all_combos,parameters,i)
    out_path = directory+names[2].format(i)
    write(result,out_path)
    return out_path
    
def multi(chunked,parameters,pool,directory,names):
    
     """Handles the multiprocessing, returning the writing paths of each individual process"""
    
    result_objs = []
    for i,chunk in enumerate(chunked):
        result=pool.apply_async(multi_numba_compiler, args=((chunk,parameters,i,directory,names)))
        result_objs.append(result)
    
    pool.close()
    pool.join()
    
    print(f"\nAlignement matrix generation ended on: {datetime.datetime.now().strftime('%c')}")
    print("\nCompiling and writing to file\n")
    
    out_paths = [result.get() for result in result_objs]
    return out_paths


#split dictionary to allocate for multiprocessing:
def dict_split(input_dict, chunks):
    
     """Divides the sgRNA dictionary into equally divided chunks corresponding
     to the total number of processes. This way each process runs a similar 
     ammount of work, only needing to be initialized once"""
    
    divider=len(input_dict)//chunks
    
    return_list = [dict() for i in range(chunks)]
    i,list_iter=0,0
    for k in input_dict:
        if i<chunks:
            return_list[i][k]=input_dict[k]
        
        if i == chunks: #odd number split will be distributed equally 
            list_iter+=1
            return_list[list_iter][k]=input_dict[k]
                
        elif len(return_list[i]) >= divider:
            i+=1

    return return_list

def write(result,directory):
    
     """writes the output of each process in the appropriate format"""
    
    with open(directory, 'w') as f:
        for guide in result:
            for matched, matrix in result[guide].items():
                f.write(guide + "\t")
                f.write(str(matrix[0]))
                for bp in matrix[1:]:
                    f.write("," + str(bp))
                f.write("\t" + matched)
                f.write("\n")
                
def write_compiled(paths,out):
    
     """Compiles all the created outputd files from the individual processes
     into one final one"""
    
    master = []
    for path in paths:
        with open(path, 'r') as f:
            for entry in f:
                master.append(entry)
        os.remove(path)
    
    with open(out, 'w') as f:
        for entry in master: 
            f.write(entry)

def parameter_parser():
    
    """ Parses the inputed length/mismatch parameters to the correct format
        input example in cmd line:
        
        7,1;9,2;20,11 
        
        organize from smallest
        to largest """
    
    parame = sys.argv[2]
    parame = parame.split(";")
    parame = [n.split(",") for n in parame]
    parame = np.array(([n for n in parame])).astype(int)
    return parame
    
def main():
    
    directory = sys.argv[1] 
    
    names = ["sg_candidates_within_genes.csv",
             "allsgcandidates.csv",
             "missmatch_matrix_{}.txt"]
    
    print(f"\nStarting at: {datetime.datetime.now().strftime('%c')}")
    
    #parameters = np.array(([7,1],[9,2],[20,11])) 
    parameters = parameter_parser()
    
    sgrna = file_parser(directory+names[0],False)
    processing_decision(sgrna,parameters,directory,names)
    
    print(f"Ended on: {datetime.datetime.now().strftime('%c')}")
    
##############
    
if __name__ == "__main__":
    multiprocessing.freeze_support() 
    main()