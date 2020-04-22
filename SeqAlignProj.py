# Sequence Alignment Problem
# Brute force vs Smith-Waterman style dynamic optimization algorithm
# Adam Cankaya
# COT 6405 Project
# Spring 2020

length_of_strings = 5   # 2 to 10,000 or more; the number of characters for each input string
number_of_runs = 3      # 1,2, or 3; the number of iterations you want to run
algorithm_number = 2    # 1 or 2; for algorithm 1 (brute force) or 2 for alogirhtm 2 (dynamic optimization)
verbose_print = True    # True or False; prints out the strings and matrices...set to False for large input string values

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
import random
from itertools import product 
import time
import numpy as np

nucleotides = ['A','C','G','T','-']
gap_penalty = 1            # point penalty for a gap
mismatch_penalty = 3       # point penalty for a mismatch
base_score = 0
optimal_score = 0
optimal_solution = list('')

# Random creates string one of length N and then
# randomly mutates it 50% to create string two
def setup(N):
    global algorithm_number, verbose_print
    str1 = ''.join(random.choice(nucleotides[0:len(nucleotides)-1]) for _ in range(N))
    str2 = ''.join(str1)

    # randomly mutate 50 percent of the characters
    mutate_rate = 0.5
    str1_list = list(str1)
    str2_list = list(str2)
    for n in range(int(N*mutate_rate)):
        random_int = random.randint(0,N-1)
        nucleotides_copy = list(nucleotides[0:len(nucleotides)-1])
        nucleotides_copy.remove(str2_list[random_int])  # don't mutate to the existing nucleotide
        str2_list[random_int] = random.choice(nucleotides_copy)
    str2 = ''.join(str2_list)

    base_score = calculate_score_alg2(str1_list, str2_list)
    if verbose_print == True:
        print("Random genome one:", str1)
        print("Random genome two:", str2)
        print('base score:', base_score)

    if algorithm_number == 1:
        return alg1(str1_list, str2_list)
    elif algorithm_number == 2:      # Call algorithm 2 - dynamic optimization 
        return alg2(str1_list, str2_list)
    else:
        print 'Something is wrong with the algorithm number'
        return 0.0


# Compares best_solution to both strings and returns the highest score
def calculate_score_alg1(best_solution, str1_list, str2_list):
    global gap_penalty, mismatch_penalty
    best_score_one = 0      # best score comparing solution to string one
    best_score_two = 0      # best score comparing solution to string two

    # if best solution is blank just compare the two strings
    if best_solution == '':
        score = 0
        for n in range(0,len(str1_list)):
            if str1_list[n] == '-' or str2_list[n] == '-':
                score += gap_penalty
            elif str1_list[n] != str2_list[n]:
                score += mismatch_penalty
        return score

    # Calculate score when best solution compared to string one
    score = 0
    for n in range(0,len(str1_list)):
        if best_solution[n] == '-':
            score += gap_penalty
        elif str1_list[n] != best_solution[n]:
            score += mismatch_penalty
    best_score_one = score

    # Calculate score when best solution compared to string two
    score = 0
    for n in range(0,len(str2_list)):
        if best_solution[n] == '-':
            score += gap_penalty
        elif str2_list[n] != best_solution[n]:
            score += mismatch_penalty
    best_score_two = score
    
    return best_score_one+best_score_two

# Compares two strings and returns the score
def calculate_score_alg2(str1_list, str2_list):
    global gap_penalty
    global mismatch_penalty
    score = 0
    for n in range(0,len(str1_list)):
        if str1_list[n] == '-' or str2_list[n] == '-':
            score += gap_penalty
        elif str1_list[n] != str2_list[n]:
            score += mismatch_penalty
    return score


# Brute force sequence alignment 
# Creates every permutation of a string of length N,
# compares them each to both string one and string two,
# prints the best results and returns the runtime
def alg1(str1_list, str2_list):
    print('Starting algorithm 1, brute force best solution')
    start_time = time.time()
    best_score = 100000
    best_solution = list('')

    # Generate list of all permutations
    for sequence in product(nucleotides,repeat=len(str1_list)):
        #print('sequence:', sequence)
        # skip sequences identical to input strings
        if ''.join(sequence) == ''.join(str1_list) or ''.join(sequence) == ''.join(str2_list):
            continue
        score = calculate_score_alg1(sequence, str1_list, str2_list)
        if score < best_score:
            best_score = score
            best_solution = list(sequence)

    end_time = time.time()
    time_elapsed = end_time - start_time    
    print('Final best score of', best_score, 'for', ''.join(best_solution))
    print('Algorithm 1 brute force completed in',"%.2f" % time_elapsed,'s')
    return time_elapsed


# Dynamic optimization sequence alignment 
# Creates every permutation of a string of length N,
# compares them each to both string one and string two,
# prints the best results and returns the runtime
def alg2(str1_list, str2_list):
    print('Starting algorithm 2, Dynamic Optimization run', runNum)
    start_time = time.time()
    score = 0
    numCols = len(str1_list)+1
    numRows = len(str2_list)+1

    # Construct scoring matrix A
    A = np.zeros((numRows,numCols))
        
    # Construct a matrix for keeping track of the arrows
    arrows = np.empty([numRows,numCols],np.dtype(object))

    # Fill bottom row left to right
    for j in range (1,numCols):
        A[numRows-1,j] = int(j * gap_penalty)
        arrows[numRows-1,j] = "left"
    
    # Fill left column bottom to top
    for i in range (numRows,0,-1):
        A[i-1,0] = int((numRows-i) * gap_penalty)
        arrows[i-1,0] = "down"
    
    # Fill scoring matrix, column by column, from the bottom to up, left to right
    # A[i,j] = min of three equations:
    #   Eq 1: A[i+1,j-1] + mismatch_penalt
    #   Eq 2: A[i,j-1] + gap_penalty
    #   Eq 3: A[i+1,j] + gap_penalty
    for j in range (1,numCols):             # for each column second from left to right
        for i in range (numRows-2,-1,-1):   # for each row second from bottom to top
            str1_index = numRows-i-2
            str2_index = j-1
            if str1_list[str1_index] != str2_list[str2_index]:    # mismatch between two strings
                eq1 = int(A[i+1,j-1] + mismatch_penalty)
            else:   # same character in both strings so no penalty
                eq1 = int(A[i+1,j-1])
            eq2 = int(A[i,j-1] + gap_penalty)
            eq3 = int(A[i+1,j] + gap_penalty)
            A[i,j] = min([eq1,eq2,eq3])
            if A[i,j] == eq1:
                arrows[i,j] = 'diag'
            elif A[i,j] == eq2:
                arrows[i,j] = 'left'
            elif A[i,j] == eq3:
                arrows[i,j] = 'down'
            else:
                print("Something broke on",i,",",j)

    if verbose_print == True:
        print(A)
        print(arrows)

    # Traceback starting at top right element of H until we hit a cell with score 0
    optimal_score = int(A[0,numCols-1])
    optimal_one = []
    optimal_two = []

    x_index = 0
    y_index = numCols-1
    count = 0
    while bool(True):
        if x_index == numRows-1 and y_index == 0:   # we're done
            break
        else:
            if arrows[x_index,y_index] == 'diag':   # both strings have same char
                str1_index = len(str1_list)-x_index-1
                str2_index = y_index-1
                optimal_one.insert(0,str1_list[str1_index])
                optimal_two.insert(0,str2_list[str2_index])
                x_index = x_index+1
                y_index = y_index-1
            elif arrows[x_index,y_index] == 'left':  # first string has a gap
                str2_index = y_index-1
                optimal_one.insert(0,'-')
                optimal_two.insert(0,str2_list[str2_index])
                y_index = y_index-1
            elif arrows[x_index,y_index] == 'down':  # second string has a gap
                str1_index = len(str1_list)-x_index-1
                optimal_one.insert(0,str1_list[str1_index])
                optimal_two.insert(0,'-')
                x_index = x_index+1
            else:   # we're done
                break
    
    end_time = time.time()
    time_elapsed = end_time - start_time
    
    print('Final optimal score of', optimal_score)
    if verbose_print == True:
        print('Final optimal alignemnts:')
        print('optimal_one:',optimal_one)
        print('optimal_two:',optimal_two)
    
    return time_elapsed


print 'Starting', number_of_runs, 'iterations of algorithm', algorithm_number, 'with input strings of length', length_of_strings
running_time = 0.0
running_total_time = 0.0
average_run_time = 0.0  
for runNum in range(1,number_of_runs+1):
    print 'Run', runNum, 'of', number_of_runs, 'for algorithm', algorithm_number, 'start'
    running_time = setup(length_of_strings)
    print 'Run', runNum, 'of', number_of_runs, 'for algorithm', algorithm_number, 'completed in', "%.2f" % running_time, 's'
    running_total_time += running_time
average_run_time = running_total_time / number_of_runs
print 'Completed', number_of_runs, 'iterations of algorithm', algorithm_number, 'with input strings of length', length_of_strings, 'in an average runtime of', "%.2f" % average_run_time, 's'
