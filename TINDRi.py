import numpy as np
from numba import njit
import multiprocessing
import datetime
from numba import types
from numba.typed import Dict
import os
import sys

# Type-expressions are currently not supported inside jit functions. 
# Leaving as global
int_array_8 = types.int8[:]
inner_dict_type = types.DictType(types.unicode_type, types.int8[:])

@njit
def numba_dict(tipo):
    
    """ Creates a numba type dictionary """
    
    dictionary = Dict.empty(
        key_type=types.unicode_type,
        value_type=tipo)
    return dictionary

@njit
def binary_subtract(array1,array2,parameters,l):
    sub = np.ones((l), dtype=np.int8)
    start_place=0
    miss=0
    for length, mismatch in parameters:
        for i,(arr1,arr2) in enumerate(zip(array1[start_place:length],array2[start_place:length])):
            if arr1-arr2 == 0:
                sub[i+start_place] = 0
            else:
                miss += 1
            if miss>mismatch:
                return #if any of the rules are broken, ignore the entire thing
        start_place=length
    return(sub)

@njit
def sgrna_all_vs_all_comparator_inner(guide,all_combos,parameters,l):
    
    """ Runs the inner loop of the all vs all sgRNA comparison.
    Compares one sgRNA to all the PAM sites. 
    Returns the mismatch matrix"""
    
    result_dict = numba_dict(int_array_8)
    for b in all_combos: 
        result = binary_subtract(guide,all_combos[b],parameters,l) 
        if result is not None:
            result_dict[b]=result

    return result_dict


def sgrna_all_vs_all_comparator_outer(chunk,all_combos,parameters,i,out_path):

    """ Runs the outer loop of the all vs all sgRNA comparison.
    Sends individually the sgRNAs for the inner loop.
    Returns the final mismatch matrix"""
    
    write_file_creator(out_path)
    for a in chunk:
        l = len(chunk[a])
        compiled = sgrna_all_vs_all_comparator_inner(chunk[a],all_combos,parameters,l)
        write(a,compiled,out_path)

def file_parser(file,nb_dict):
    
    """ Parses the input .csv files into dictionaries. Converts all DNA
    sequences to their respective binary array forms. This gives some computing
    speed advantages."""
    
    if nb_dict:
        container = numba_dict(int_array_8)
    else:
        container = {}
    
    with open(file) as current: 
        for line in current:
            line = line[:-1].split(",")
            name = line[0]
            sequence = line[1].upper()
            sequence = sequence.replace(" ", "")
            byte_list = bytearray(sequence,'utf8')
            if name not in container:
                container[name] = np.array((byte_list), dtype=np.int8)
                
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
    chunked = dict_split(sgrna, cpu)
    out_paths=multi(chunked,parameters,pool,directory,names)
    out = directory+names[2].format("end")
    write_compiled(out_paths,out)

def multi_numba_compiler(chunk,parameters,i,directory,names):
    
    """Initializes the numba dictionaries with all the sgRNAs, and PAM sites.
    Because numba objects cannot be pickled at the moment, 
    all numba processing needs to be done after child process is created"""
    
    all_combos = file_parser(directory+names[1],True)
    chunk_numba = numba_dict(int_array_8)
    out_path = directory+names[2].format(i)

    for chunky in chunk:
        chunk_numba[chunky] = chunk[chunky]
    
    sgrna_all_vs_all_comparator_outer(chunk_numba,all_combos,parameters,i,out_path)
    return out_path
    
def multi(chunked,parameters,pool,directory,names):
    
    """Handles the multiprocessing, 
     returning the writing paths of each individual process """
    
    result_objs = []
    for i,chunk in enumerate(chunked):
        result=pool.apply_async(multi_numba_compiler, args=((chunk,parameters,i,directory,names)))
        result_objs.append(result)
    
    pool.close()
    pool.join()
    
    print(f"\nAlignment matrix generation ended on: {datetime.datetime.now().strftime('%c')}")

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

def write_file_creator(out_path):
    if not os.path.isfile(out_path):
        with open(out_path, "w") as text_file:
            text_file.close()

def write(guide,result,out_path):
    
    """writes the output of each process in the appropriate format"""

    with open(out_path, "a+") as f:
        for matched in result:
            matrix = result[matched]
            f.write(guide + "\t")
            f.write(str(matrix[0]))
            for bp in matrix[1:]:
                f.write("," + str(bp))
            f.write("\t" + matched)
            f.write("\n")

def write_compiled(paths,out):
    
    """Compiles all the created outputd files from the individual processes
     into one final one"""
    write_file_creator(out)

    master = []
    for path in paths:
        with open(path, 'r') as f:
            for entry in f:
                master.append(entry)
        os.remove(path)
        with open(out, 'a+') as f:
            for entry in master: 
                f.write(entry)
        master = []

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
    
    directory = sys.argv[1] #just the directory, as the names are hardcoded

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
