from pylab import *
import pickle
from os.path import isfile
from scipy.linalg import expm

# Time to collect K coupons of
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

def ODE_K_coupon(n,k):
    x  = array([0.]*(k+1)); x[0] = 1;
    dx = array([0.]*(k+1));
    h = 0.00001;
    traj = [0]*int(20/h);
    for t in range(0,int(20/h)):
        dx[0] = 0;
        for i in range(0,k):
            dx[i] -= h*x[i];
            dx[i+1] = h*x[i];
        x = x+dx;
        traj[t] = x[k];
    t = array(arange(0,20,h))*n;
    return(t,traj)

def t_N_K_coupon_ODE(n,k):
    (t,traj) = ODE_K_coupon(n,k)
    return(t[list(array(traj)>1-1./n).index(True)])

def t_N_K_coupon_exp(n,k):
    Q = matrix( [[0.]*k for i in range(0,k)] );
    for i in range(0,k):
        Q[i,i] = -1;
    for i in range(0,k-1):
        Q[i+1,i] = 1;
    t = 5.;
    while sum(expm(Q*t)[k-1,:]) > 1./n: t *= 1.5;
    dt = t*0.8;
    for i in range(0,30):
        if sum(expm(Q*t)[k-1,:]) > 1./n:
            t+= dt;
        else:
            t -= dt;
        dt /= 2;
    return (n*t)

def t_N_K_coupon(n,k):
    return(t_N_K_coupon_exp(n,k))

gamma_constant  = 0.5770836328

n = 100;
k = 5;
nb_samples = 10000

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


def plot_figures(n,k):
    (t_ODE,traj_ODE) = ODE_K_coupon(n,k);
    t_N = t_N_K_coupon(n,k);
    
    f = figure(1);
    clf();
    plot(t_ODE,traj_ODE,'r-')
    for i in range(0,10):
        traj_sto = traj_K_coupon(n,k);
        plot(traj_sto,'k+--')
    plot([t_N,t_N],[0,1])
    text(t_N+5,0.1,'t_N')
    f.savefig('traj_compa_n{}_k{}.pdf'.format(n,k));
        
        
    f=figure(2);
    clf();
    dist = compute_or_load_dist(n,k,nb_samples)
    #dist = [len(traj_K_coupon(n,k)) for i in range(0,nb_samples)]
    a = histogram(dist,50)
    plot(a[1][1:],cumsum(a[0])/nb_samples)
    plot([t_N + n*x for x in arange(-5,5,.1)],
         [exp(-exp(-x)) for x in arange(-5,5,.1)],'--')
    plot([t_N + n*x*log(k) for x in arange(-5,5,.1)],
         [exp(-exp(-x)) for x in arange(-5,5,.1)],'--')
    legend(('Simulation ($n={}$)'.format(n),'Gumbel approximation','Gumbel approximation (log(k))'),loc='best')
    print('n={},k={}, Mean of T_N = {:.2f}+/-{:.2f}, t_N+n*\gamma={:.2f}, *log(k)={:.2f}'.format(
        n,k,mean(dist), 2*sqrt(var(dist)/nb_samples),
        t_N+gamma_constant*n, t_N+gamma_constant*n*log(k)))
    f.savefig('hist_compa_n{}_k{}.pdf'.format(n,k));

    
# for n in [10,20,50,100]:
#     for k in [1,2,5,10]:
#         plot_figures(n,k);

# figure(3);
# clf();
# dist = [len(traj_K_coupon(n,k)) for i in range(0,10000)]
# a = histogram(dist,50)
# plot(a[1][1:],cumsum(a[0]),)
# plot([t_N + n*x for x in arange(-5,5,.1)],
#      [exp(-exp(-x)) for x in arange(-5,5,.1)],'--')


# mean([len(traj_K_coupon(n,k)) for i in range(0,100000)])

ls = ['-', ':', '-.' , ':' , '--']

for k in [1,5]:
    f = figure(1);
    f.set_size_inches(6,4.5)
    clf();
    if k==1:
        my_n = [10,100];
    else:
        my_n = [10,100,10000];
    for (i,n) in enumerate(my_n):
        dist = compute_or_load_dist(n,k)
        t_N = t_N_K_coupon(n,k);
        plot( (sort(dist)-t_N)/n, arange(0,1,0.0001), ls[i])
        print('n={},k={}, \tMean of T_N = {:.2f}\t+/-{:.2f},\tt_N+n*\gamma={:.2f}'.format(
            n,k,mean(dist), 2*sqrt(var(dist)/nb_samples),
            t_N+gamma_constant*n))#, t_N+gamma_constant*n*log(k)))
    #plot((a[1][1:]-t_N)/n,cumsum(a[0])/nb_samples)
    plot([x for x in arange(-5,5,.1)],
         [exp(-exp(-x)) for x in arange(-5,5,.1)],'--')
    xlim(-2,4)
    xlabel('$(T_n-t_n)/n$')
    ylabel('cumulative probability')
    if k==1:
        legend((r'$n=10$',r'$n=100$','Gumbel distribution'),loc='best')
    else:
        legend((r'$n=10$',r'$n=100$',r'$n=10000$','Gumbel dist.'),loc='best')
    f.savefig('gumbel_k{}.pdf'.format(k),bbox_inches='tight');

    
# for k in [1,2,5,10]:
#     for n in [10,100,500,1000,5000,10000]:
#         dist = compute_or_load_dist(n,k)
#         print('{}-n={}   \t{}'.format(k,n,(mean(dist)-t_N_K_coupon(n,k))/n-gamma_constant))
        
