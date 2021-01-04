import numpy as np
import math


def sum_up(matrix, row):
    res=0
    for i in range(len(matrix)):
        res=res+matrix[row][i]
        
    return res

def compare(vector, number):
    res=0
    for x in vector:
        if number<=vector[res]:
            return res
        else:
            res=res+1
    return -1

def rate_interval(matrix, row):
    vec=[]
    for col in range(1,len(matrix)+1):
        sum=0
        for i in range(col):
            sum=sum+matrix[row][i]
        res=sum/sum_up(matrix, row)
        vec.append(res)
    #print (vec)
    return vec


#initialisierung 
def init(matrix):
    sumUp=[]
    rateInterval=[]
    for row in range(len(matrix)):
        sumUp.append(sum_up(matrix,row))
        rateInterval.append(rate_interval(matrix,row))
    return [sumUp, rateInterval]
    #print(sumUp)
    #print(rateInterval)



def timestep(oldstate, t1, sumUp, rateInterval):
    dt=-np.log(np.random.random_sample() )/sumUp[oldstate]
    #print(dt)
    t1=t1+dt
    z1=np.random.random_sample()
    #print(z1)
    newstate=compare(rateInterval[oldstate], z1)
    #print(newstate)
    return [newstate, t1]

def simulation_prob(matrix, startState, T):
    [sumUp, rateInterval]=init(matrix)
    t=0.0
    state=startState
    #memory for trajectory
    TIME=[]
    states=[]

    while t<T:
        [state, t]=timestep(state, t, sumUp, rateInterval)
        TIME.append(t)
        states.append(state)
    prob=np.zeros(len(matrix))

    for i in range(len(TIME)-1):
        t=TIME[i+1]-TIME[i]
        prob[ states[i] ]=prob[ states[i] ]+t
    prob=1.0/T*prob
    return prob

def simulation_entropy_prod(matrix, startState, T):
    #start_time = time.time()
    [sumUp, rateInterval]=init(matrix)
    t=0.0
    state=startState
    W=0
    while t<T:
        oldstate=state
        [state, t]=timestep(state, t, sumUp, rateInterval)
        W=W+math.log(matrix[oldstate][state]/matrix[state][oldstate])
    #print("--- %s seconds ---" % (time.time() - start_time))
    #print W/T
    return W/T

def pdf_entropy(Matrix, start, time, runs):
    entropy=[]
    for i in range(runs):
        W=simulation_entropy_prod(matrix=Matrix, startState=start, T=time)
        entropy.append(W)
    return entropy

def pdf_flux(Matrix, start, time, runs, state1, state2):
    flux=[]
    for i in range(runs):
        p=simulation_prob(matrix=Matrix, startState=start, T=time)
        J=p[state1]*Matrix[state1][state2]-p[state2]*Matrix[state2][state1]
        flux.append(J)
    return flux

def simulation_transitions(matrix, startState, T, state1, state2):
    total_transitions=0
    transition_12=0
    transition_21=0
    [sumUp, rateInterval]=init(matrix)
    t=0.0
    state=startState
    while t<T:
        oldstate=state
        [state, t]=timestep(state, t, sumUp, rateInterval)
        total_transitions=total_transitions+1
        if oldstate==state1 and state==state2:
            transition_12=transition_12+1
        elif oldstate==state2 and state==state1:
            transition_21=transition_21+1
    return [transition_12, transition_21, total_transitions]

def pdf_transition(Matrix, start, time, runs, state1, state2):
    totalcounts=[]
    counts12=[]
    counts21=[]
    netcounts=[]
    for i in range(runs):
        [t12, t21, total]=simulation_transitions(matrix=Matrix, startState=start, T=time, state1=state1, state2=state2)
        totalcounts.append(total)
        counts12.append(float(t12))
        counts21.append(float(t21))
        netcounts.append(float(t12-t21))
    totalcounts=np.array(totalcounts)*1.0/time
    counts12=np.array(counts12)*1.0/time
    counts21=np.array(counts21)*1.0/time
    netcounts=np.array(netcounts)*1.0/time
    return [totalcounts, counts12, counts21, netcounts]
    #return netcounts

def simulation_both(matrix, startState, T, state1, state2):
    total_transitions=0
    transition_12=0
    transition_21=0
    [sumUp, rateInterval]=init(matrix)
    t=0.0
    state=startState
    W=0
    while t<T:
        oldstate=state
        [state, t]=timestep(state, t, sumUp, rateInterval)
        W=W+math.log(matrix[oldstate][state]/matrix[state][oldstate])
        total_transitions=total_transitions+1
        if oldstate==state1 and state==state2:
            transition_12=transition_12+1
        elif oldstate==state2 and state==state1:
            transition_21=transition_21+1
    return transition_12, transition_21, total_transitions, W/T   
    
def pdf_transitions_and_entropy(Matrix, start, time, runs, state1, state2):
    #totalcounts=[]
    #counts12=[]
    #counts21=[]
    netcounts=[]
    
    entropy=[]

    for i in range(runs):
        t12, t21, total, P=simulation_both(matrix=Matrix, startState=start, T=time, state1=state1, state2=state2)
        entropy.append(P)
        #totalcounts.append(total)
        #counts12.append(float(t12))
        #counts21.append(float(t21))
        netcounts.append(float(t12-t21))
    #totalcounts=np.array(totalcounts)*1.0/time
    #counts12=np.array(counts12)*1.0/time
    #counts21=np.array(counts21)*1.0/time
    netcounts=np.array(netcounts)*1.0/time
    return netcounts, entropy
    
