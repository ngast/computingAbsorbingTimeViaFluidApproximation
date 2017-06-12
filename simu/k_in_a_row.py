from pylab import *
from os.path import isfile
import pickle
from scipy.linalg import expm


gamma_constant  = 0.5770836328

def m_approx(Q,t,proba_dist,n):
    if proba_dist is None:
        return(sum(expm(Q*t/n)[-1,:]))
    else:
        return(mean([sum(expm(Q*t*p)[-1,:]) for p in proba_dist]))
def compute_or_load_dist(n,k, new=False, nb_samples=10000):
    if (isfile('dict_dist_N_K_coupon')):
        with open('dict_dist_N_K_coupon','rb') as f: dict_coupon = pickle.load(f);
    else:
        dict_coupon = dict([])
    if (not ((n,k) in dict_coupon)) or new:
        print('computing dist for n={} k={}'.format(n,k));
        dist = [len(traj_K_coupon(n,k)) for i in range(0,nb_samples)]
        if (isfile('dict_dist_N_K_coupon')):
            with open('dict_dist_N_K_coupon','rb') as f: dict_coupon = pickle.load(f);
        else:
            dict_coupon = dict([])

        dict_coupon[(n,k)] = dist
        with open('dict_dist_N_K_coupon','wb') as f:
            pickle.dump(dict_coupon,f,pickle.HIGHEST_PROTOCOL)
    return(dict_coupon[(n,k)])
def traj_K_coupon(n,k=1):
    coupons = [0]*n;
    traj = [];
    nb_finished = 0;
    while(nb_finished < n): 
        j = randint(n); 
        coupons[j] += 1
        if coupons[j] == k:
            nb_finished += 1;
        traj.append(nb_finished/n)
    return(traj)

    
def compute_or_load_dist_k_in_a_row(n,k,p, new=False, proba_dist=None,nb_samples=10000):
    if (isfile('dict_dist_N_K_in_a_row')):
        with open('dict_dist_N_K_in_a_row','rb') as f: dict_in_a_row = pickle.load(f);
    else:
        dict_in_a_row = dict([])
    if proba_dist == None:
        my_name = (n,k,p)
    else:
        my_name = (n,k,p,proba_dist,nb_samples)
    if (not (my_name in dict_in_a_row)) or new:
        print('computing dist for n={} k={}'.format(n,k));
        print(my_name)
        dist = [len(traj_K_in_a_row(n,k,p,proba_dist=proba_dist)) for i in range(0,nb_samples)]
        if (isfile('dict_dist_N_K_in_a_row')):
            with open('dict_dist_N_K_in_a_row','rb') as f: dict_in_a_row = pickle.load(f);
        else:
            dict_in_a_row = dict([])
        dict_in_a_row[my_name] = dist
        with open('dict_dist_N_K_in_a_row','wb') as f:
            pickle.dump(dict_in_a_row,f,pickle.HIGHEST_PROTOCOL)
    return(dict_in_a_row[my_name])

# Time to collect K in_a_rows of
def traj_K_in_a_row(n,k,p,proba_dist=None):
    if proba_dist == 'log':
        proba_dist = array([1/log(1+i) for i in range(1,n+1)])
        proba_dist = proba_dist / sum(proba_dist);
    elif proba_dist == 'square':
        proba_dist = array([1/(i*i) for i in range(1,n+1)])
        proba_dist = proba_dist / sum(proba_dist);
    elif proba_dist != None:
        print('other distribs non implemented')
    coupons = [0]*n;
    traj = [];
    nb_finished = 0;
    while(nb_finished < n): 
        j = choice(n,p=proba_dist);
        if coupons[j] < k:
            if rand() < p:
                coupons[j] += 1
            else:
                coupons[j] = 0;
            if coupons[j] == k:
                nb_finished += 1;
        traj.append(nb_finished/n)
    return(traj)

def Q_for_K_in_a_row(n,k,p):
    Q = matrix( [[0.]*k for i in range(0,k)] );
    for i in range(0,k):
        Q[i,i] = -1;
        Q[i,-1] += 1-p
    for i in range(0,k-1):
        Q[i+1,i] = p;
    return(Q)
def t_N_K_in_a_row(n,k,p,proba_dist=None):
    if proba_dist == 'log':
        proba_dist = array([1/log(1+i) for i in range(1,n+1)])
        proba_dist = proba_dist / sum(proba_dist);
    elif proba_dist == 'square':
        proba_dist = array([1/(i*i) for i in range(1,n+1)])
        proba_dist = proba_dist / sum(proba_dist);
    elif proba_dist!=None:
        print('other distribs non implemented')
    Q = Q_for_K_in_a_row(n,k,p)
    t = 5*n;
    #print(m_approx(Q,t,proba_dist),(Q,t,proba_dist))
    while m_approx(Q,t,proba_dist,n) > 1./n: t *= 1.5;
    dt = 0.5*t;
    for i in range(0,30):
        #print(t)
        if m_approx(Q,t,proba_dist,n) > 1./n:
            t+= dt;
        else:
            t -= dt;
        dt /= 2;
    return (t)

def T_N_K_in_a_row_Approx(n,k,p,proba_dist=None):
    t_N =  t_N_K_in_a_row(n,k,p,proba_dist=proba_dist)
    Q = Q_for_K_in_a_row(n,k,p)
    multiply = 1/min(abs(eig(Q)[0]));
    #multiply = -1/max(eig(Q)[0]);
    #print(t_N,multiply)
    if proba_dist=='log': multiply *= sum([1/log(1+i) for i in range(1,n+1)])*log(n)
    elif proba_dist=='square': multiply *= sum([1/(i*i) for i in range(1,n+1)])*(n+1)*(n+1)
    elif proba_dist ==None: multiply *= n;
    return(t_N + multiply*gamma_constant)

def borne(n,k,p):
    return( n*(log(n)+gamma_constant)*(p**(-k)-1)/(1-p) )

def simu(n,k,p):
    coupon = 0;
    t=0;
    while( coupon < k ):
        if rand() < p:
            coupon += 1
        else:
            coupon = 0;
        t += 1
    return(t);

def perf_function_p(n,k,my_p):
    bor = [0]*len(my_p);
    t_N = [0]*len(my_p);
    T_N = [0]*len(my_p);
    conf_T_N = [0]*len(my_p);
    for i,p in enumerate(my_p):
        dist = compute_or_load_dist_k_in_a_row(n,k,p)
        bor[i] = borne(n,k,p)
        t_N[i] = real(T_N_K_in_a_row_Approx(n,k,p))
        T_N[i] = mean(dist)
        conf_T_N[i] = sqrt(var(dist)/len(dist))
        print(p,borne(n,k,p), real(T_N_K_in_a_row_Approx(n,k,p)),mean(dist))
    return(bor,t_N,T_N,conf_T_N)


k = 10
n = 10
my_p = arange(0.70,1,.025)

f = figure(1);
f.set_size_inches(6,4.5)
clf()
for n in [10]:
    bor,t_N,T_N,conf_T_N = perf_function_p(n,k,my_p)
    plot(my_p,bor,'k--')
    plot(my_p,t_N,'kx-')
    errorbar(my_p,T_N,conf_T_N)
legend((r'Bound $(T_1^\max n(\log n+\gamma)$)',r'Approximation ($t_n+n\gamma/\nu$)','Simulation'))
xlabel('Parameter $q$')
ylabel('Hitting time')
#text(0.705, 2000, r'$n=10$')
#text(0.82, 2700, r'$n=20$')
xlim(0.7,1)
ylim(0,3200)
f.savefig('erasure_channel.pdf',bbox_inches='tight');


# p = .5
    
# for k in [1,2,5,10]:
#     for n in [10,20,50,100]:
#         dist = compute_or_load_dist_k_in_a_row(n,k,p)
#         dist_c = compute_or_load_dist(n,k)
#         Q = Q_for_K_in_a_row(n,k,p);
#         print(mean(dist), mean(dist_c), t_N_K_in_a_row(n,k,p)+n*gamma_constant*(-1/max(eig(Q)[0])), T_N_K_in_a_row_Approx(n,k,p))
#         print('{}-n={}   \t{}'.format(k,n,(mean(dist)-t_N_K_in_a_row(n,k,p))/n-gamma_constant))

