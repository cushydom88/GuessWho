### Probabilities from recursion formula ###

# For board state (N, M) compute the probability of winning with optimal play

# When computing a large number of probabilities call gw_matrix once and pass the result in to gw_prob_using_matrix multiple times
# Compute the matrix (a_i,j) where a_i,j = i*j*p(i, j)
def gw_matrix( maxCards = 32, Ternary = False ):
    A_arr = [ [  -1  for m in range( maxCards) ] for n in range( maxCards ) ]
    for m in range(1, maxCards+1) :
        A_arr[0][m-1] = m
    A_arr[1][0] = 1
    A_arr[1][1] = 2
    for n in range(3, maxCards+1):
        for m in range(1, n):
            if Ternary:
                my_arr = [ (m*n - A_arr[m-1][t-1] -  A_arr[m-1][n-t-1] ) for t in range(1,1+(n//2) )  ]
                for s in range(1,1+(n//3) ):
                    my_arr.extend( [ (m*n - A_arr[m-1][s-1] - A_arr[m-1][t-1] -  A_arr[m-1][n-s-t-1] ) for t in range(s,1+((n-s)//2) )  ]  )
            else:
                my_arr = [ (m*n - A_arr[m-1][s-1] -  A_arr[m-1][n-s-1] ) for s in range(1,1+(n//2))  ] 
            my_arr.append(m)
            MyMax = max(my_arr)
            A_arr[n-1][m-1] = MyMax
        for k in range(2,n+1):
            if Ternary:
                my_arr = [(k*n - A_arr[n-1][t-1] -  A_arr[n-1][k-t-1] ) for t in range(1,1+(k//2) ) ]
                for s in range(1,1+(k//3) ):
                    my_arr.extend( [ (k*n - A_arr[n-1][s-1] - A_arr[n-1][t-1] -  A_arr[n-1][k-s-t-1] ) for t in range(s,1+((k-s)//2) )  ]  )
            else:
                my_arr = [ (k*n - A_arr[n-1][s-1] -  A_arr[n-1][k-s-1] ) for s in range(1,1+(k//2))  ] 
            my_arr.append(n)
            MyMax = max(my_arr)
            A_arr[k-1][n-1] = MyMax
    return A_arr

def gw_prob_using_matrix(N, M, Mat ):
    maxCards = len(Mat)
    if (N>maxCards) or (M>maxCards):
        print('Error: Invalid values for N, M.')
        return 
    return (100.0*Mat[N-1][M-1])/(1.0*N*M)

# Call the following function when processing a small number of board states
def gw_prob(N, M, Ternary = False):
    maxCards = max([N, M])
    Mat = gw_matrix( maxCards, Ternary )
    return gw_prob_using_matrix(N, M, Mat )


### Optimal Strat Functions

def find_p(n):
    P = 1
    while n > 3**P:
        P +=1
    return P

# Define an optimal strategy
def Strat(A, B, Ternary = False ):
    if Ternary:
        if A == 1:
            return [0,0,1]
        elif A == 2:
            return [0,1,1]
        elif A == 3:
            return [1,1,1]
        elif A == 4:
            return [1,1,2]
        elif A == 5:
            if B <= 4:
                return [1,2,2]
            else:
                return [1,1,3]
        elif A == 6:
            if B <= 4:
                return [2,2,2]
            elif B <= 7:
                return [1,1,4]
            else:
                return [1,2,3]
        elif A == 7:
            if B <= 4:
                return [2,2,3]
            else:
                return [1,3,3]
        elif A == 8:
            if B <= 4:
                return [2,3,3]
            elif B <= 7:
                return [1,3,4]
            else:
                return [2,3,3]
        elif A == 9:
            if B <= 4:
                return [3,3,3]
            elif B <= 7:
                return [1,4,4]
            else:
                return [3,3,3]
        else:
            C = find_p(A)
            E = 3**(C-1)
            E1 = 3**(C-2)
            E2 = 3**(C-3)
            if B < E + E1:
                L1 = A - 3*E1
                R = L1//(3*E2)
                X1 = E1 + R*E2
                L2 = A - 3*X1
                Z = X1 + min([L2,E2])
                L3 = A - 2*X1 - Z
                Y = X1 + min([L3,E2])
                L4 = A - X1 - Y - Z
                X = X1 + L4
                Res = [X,Y,Z]
                return Res
            else:
                if (A==(1+E + 2*E1)) and (B > 2*E) and (B < (3**C)-E1):
                    return [E1,E1+2,E-1]
                elif (A==(1+ 2*E + E1)) and (B > 2*E) and (B < (3**C)-E1):
                    return [E1+2,E-1,E]
                elif A <= E + 2*E1:
                    Res = [E1, E1, A-2*E1]
                elif A <= 2*E + E1:
                    Res = [E1, A-E-E1, E ]
                else:
                    Res = [A-2*E, E,E ]
                return Res
    else:
        if A == 1:
            return [0,1]
        elif A == 2:
            return [1,1]
        elif (A == 4) and (B == 4):
            return [1,3]
        elif (A == 6) and (B == 4):
            return [3,3]
        elif (A == 10) and (B == 4):
            return [5,5]
        elif (A%4) == 0:
            return [A//2,A//2]
        elif (A%4) == 2:
            return [(A//2)-1,(A//2)+1]
        else:
            return [A//2,(A//2)+1]

def Use_Strat_1(N1, M1):
    Res = 0
    Stack = [[N1,M1]]
    Counts = [1]
    while Stack != []:
        Curr = Stack[0]
        N, M = Curr[0], Curr[1]
        V = Counts[0]
        StackN = Stack[1:]
        CountsN = Counts[1:]
        if N == 1:
            Res += M*V
        elif M == 1:
            Res += V
        else:
            Res += N*M*V
            S = Strat(N,M)
            New1, New2 = [M, S[0]], [M, S[1]]
            if New1 in StackN:
                I = StackN.index(New1)
                CountsN[I] = CountsN[I] - V
            else:
                StackN.append(New1)
                CountsN.append(-V)
            if New2 in StackN:
                I = StackN.index(New2)
                CountsN[I] = CountsN[I] - V
            else:
                StackN.append(New2)
                CountsN.append(-V)
            StackN = [ StackN[i] for i in range(len(StackN)) if CountsN[i] !=0 ]
            CountsN = [ CountsN[i] for i in range(len(CountsN)) if CountsN[i] !=0 ]
        Stack = StackN
        Counts = CountsN
    return Res

def tern_prob_small( A, B ):    
    if B == 1:
        return 1
    else:
        ARR = [[4, 7, 10, 13, 16, 20, 24, 28], [5, 8, 11, 15, 19, 24, 29, 34], [6, 9, 12, 17, 22, 27, 32, 38], [6, 9, 13, 18, 23, 30, 37, 44], [6, 9, 14, 20, 26, 33, 40, 48], [6, 9, 15, 22, 29, 36, 45, 54]]
        return(ARR[A-4][B-2])

def Use_Strat_2(N1, M1):
    Res = 0
    Stack = [[N1,M1]]
    Counts = [1]
    while Stack != []:
        Curr = Stack[0]
        N, M = Curr[0], Curr[1]
        V = Counts[0]
        StackN = Stack[1:]
        CountsN = Counts[1:]
        if N == 1:
            Res += M*V
        elif N == 2:
            if M == 1:
                Res += V
            else:
                Res += 2*(M - 1)*V
        elif N == 3:
            if M == 1:
                Res += V
            else:
                Res += 3*(M-1)*V
        elif ( N >= 4 ) and ( N <= 9 ):
            if M <= 9:
                Val = tern_prob_small(N, M)
                Res += V*Val
            else:
                C = 2
                Val = N*M + 2*(3**(2*C - 1))  - (3**(C-2))*( 6*( (1 + 3**C - N)//2  ) + 4*( ( 3**C - N)//2  ) )- N*(3**C)
                Res += V*Val
        elif M == 1:
            Res += V
        else:
            Res += N*M*V
            S = Strat(N,M)
            New1, New2, New3 = [M, S[0]], [M, S[1]], [M,S[2]]
            if New1 in StackN:
                I = StackN.index(New1)
                CountsN[I] = CountsN[I] - V
            else:
                StackN.append(New1)
                CountsN.append(-V)
            if New2 in StackN:
                I = StackN.index(New2)
                CountsN[I] = CountsN[I] - V
            else:
                StackN.append(New2)
                CountsN.append(-V)
            if New3 in StackN:
                I = StackN.index(New3)
                CountsN[I] = CountsN[I] - V
            else:
                StackN.append(New3)
                CountsN.append(-V)
            StackN = [ StackN[i] for i in range(len(StackN)) if CountsN[i] !=0 ]
            CountsN = [ CountsN[i] for i in range(len(CountsN)) if CountsN[i] !=0 ]
        Stack = StackN
        Counts = CountsN
    return Res

# For board state (N, M) compute the probability of winning with optimal play by using optimal play 
def gw_prob_from_strat(N, M, Ternary = False):
    if Ternary:
        return (100.0*Use_Strat_2(N,M))/(1.0*N*M)
    else:
        return (100.0*Use_Strat_1(N,M))/(1.0*N*M)


# Verify the optimal strategy using the recursion relations

def verify_strat( N, M, Mat = None, Ternary = None ):
    if Ternary is None:
        Ternary = False
    if Mat is None:
        maxCards = max([N, M])
        Mat = gw_matrix( maxCards, Ternary )
    S = Strat(N, M, Ternary )
    E = N*M
    for x in S:
        E -= Mat[M-1][x-1]
    return Mat[N-1][M-1] == E


### Example ###
Pass = True
Mat = gw_matrix()
for N in range(3,33):
    for M in range(3,33):
        if not verify_strat(N, M, Mat):
            print(N, M, ' Failed')
            Pass = False
if Pass:
    print('All Passed')