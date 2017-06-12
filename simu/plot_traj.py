from pylab import *

# To regenerate the data
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

seed(7)

a = traj_K_coupon(200,3)
s200 = array([[i/200,a[i]] for i in range(len(a))])
a = traj_K_coupon(20,3)
s20 = array([[i/20,a[i]] for i in range(len(a))])

def T_n(s):
    i = 0
    while s[i,1] < 1:
        i+=1
    return i

t = arange(0,70,.01)
ex = [1-(1+x+x*x/2)*exp(-x) for x in t]
T = T_n(s20)

f = figure(1);
f.set_size_inches(6,5)
clf()
plot(20*s20[:,0], s20[:,1],'+-')
plot(20*t,ex)
plot([0,7*20],[.95,.95],'-.')
plot([125.9,125.9],[0,1],'k--')
plot([T,T],[0,1],'k--')
text(T-8,.05,r'$T_n$')
text(128,.05,r'$t_n$')
yticks([0,0.2,0.4,0.6,0.8,0.95,1],[0,0.2,0.4,0.6,0.8,r'$1-\frac{1}{n}$'])
xlim(0,7*20)
xlabel('time')
ylabel(r'$M_0(t)$')
legend(('Simulation','Fluid approximation'),loc='upper left')
f.savefig('simu_20.pdf',bbox_inches='tight');


# f = figure(1);
# f.set_size_inches(3,5)
# clf()
# plot(20*s20[:,0], s20[:,1],'+-')
# plot(20*t,ex)
# plot([0,7*20],[.95,.95],'-.')
# plot([125.9,125.9],[0,1],'k--')
# plot([121,121],[0,1],'k--')
# text(113,.55,r'$T_n$')
# text(128,.55,r'$t_n$')
# text(-14,.93,r'$1-\frac{1}{n}$')
# xlim(60,140)
# ylim(.5,1)
# xlabel('time')
# ylabel(r'$M_0(t)$')
# f.savefig('simu_20_zoom.pdf',bbox_inches='tight');



f = figure(1);
f.set_size_inches(6,5)
clf()
plot(200*s200[:,0], s200[:,1],'+-')
plot(200*t,ex)
plot([0,10*200],[.995,.995],'-.')
plot([1854.8,1854.8],[0,1],'k--')
T = T_n(s200)
plot([T,T],[0,1],'k--')
text(1775,0.06,r'$t_n$')
text(T+14,0.06,r'$T_n$')

#arrow(1750,0,10,0.1, head_width=.15, head_length=10);
#text(1862,.05,r'$t_n$')
xlabel('time')
ylabel(r'$M_0(t)$')
xlim(0,10*200)
yticks([0,0.2,0.4,0.6,0.8,0.995,1],[0,0.2,0.4,0.6,0.8,r'$1-\frac{1}{n}$'])
legend(('Simulation','Fluid approximation'),loc='best')
f.savefig('simu_200.pdf',bbox_inches='tight');

f = figure(1);
f.set_size_inches(3,5)
clf()
plot(200*s200[:,0], s200[:,1],'+-')
plot(200*t,ex)
plot([0,10*200],[.995,.995],'-.')
plot([1854.8,1854.8],[0,1],'k--')
plot([T,T],[0,1],'k--')
xlabel('time')
ylabel(r'$M_0(t)$')
xlim(1500,1950)
ylim(0.95,1)
xticks([1500,1600,1700,1800,1900])
yticks([0.95,0.96,0.97,0.98,0.99,0.995],[0.95,0.96,0.97,0.98,0.99,r'$1-\frac{1}{n}$'])
text(T+9,0.96,r'$T_n$')
text(1810,.96,r'$t_n$')
f.savefig('simu_200_zoom.pdf',bbox_inches='tight');
