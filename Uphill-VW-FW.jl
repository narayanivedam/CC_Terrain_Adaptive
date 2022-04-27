
using PyPlot

#Discretization:
T = 2
h = 50

#Engine Map, Empirical Constants:
a_e = -0.1135
b_e = 9.4263
c_e = 5.6282
a_d = -0.0917
b_d = -46.0140

#Longitudinal Dynamics Parameters:
c_w = 0.6#Drag Coefficient
r_w = 0.52#Wheel Radius
A = 10#Frontal Area
rho = 1.29#Density Of Air
mass = 40000#Vehicle Mass
g = 9.8#Acceleration due to gravity
J_w = 32.9#Wheel Inertia
k_b = 20e3#Total Brake force
J_e = 3.5#Inertia Of the Engine
n_cyl = 6#No of engine cylinders
n_r = 2 #No. of revolutions per stroke
c_r = 7e-3#Rolling Resistance Coefficient

#Fuel Flow, Empirical Constants:
a_del = -1.429e-4
b_del = 0.3973
c_del = -48.5649

#Transmission Parameters:
eta_t = 0.97#Transmission efficiency
eta_f = 0.97#Final drive efficiency
i_f = 3.27#Final drive ratio
i_g = [11.27; 9.14; 7.17; 5.81; 4.62; 3.75; 3.01; 2.44; 1.91; 1.55; 1.23; 1]#Transmission ratio as a function of gear position

#Constants:
c_1 = 0
c_2 = 0.5*rho*c_w*A*r_w#Air Drag
c_3 = 0
c_7 = atan(c_r)
c_8 = mass*g*r_w*sqrt(1+(c_r*c_r))
c_5 = n_cyl/(6e4*n_r)

#Initialization:
#State-Space Discretization
v = collect(80:0.1:90)
v_min = minimum(v)
v_max = maximum(v)
v_ref = 85
m = length(v)
k2 = collect(Int,1:1:51)
n = length(k2)
p = zeros(n,m,m)
i2 = collect(Int,1:1:m)
j2 = collect(Int,1:1:m)
b = zeros(n,m,m)
gr = zeros(n,m,m)
mf_dot = 0
m_f = zeros(n+1,m,m)
m_f[1,:,:] = 5
J = zeros(n,m,m)
del_max = zeros(m,1)#Maximum Fueling
e = zeros(m,1)
i_t = zeros(m,1)
N = zeros(m+1,1)
G = zeros(m+1,1)
#Road slope
alpha = zeros(n,1)
alpha[5:26] = 0.03
u = zeros(m,1)

#Intial Conditions
N[1] = (30*minimum(v)*i_f)/(3.6*pi*r_w)
G[1] = 12

#Cost Function Parameters:
q1=10
q2=5
q3=15
q4=15
q5=100
del = 0.1
i = 1
j = 1
k = 1

#Stage
for k in k2
    
    #State in stage k
    for i in i2
        
        #Gear Shift
        if(G[i]==1)
            if(N[i]>=1500)
                G[i+1]= 2
                else
                G[i+1] = 1
            end
            
        elseif(G[i]==2)
            if(N[i]>=1501)
                G[i+1] = 3
            elseif(N[i]<=950)
                G[i+1]= 1
            else
                G[i+1] = 2
            end
            
        elseif(G[i]==3)
            if(N[i]>=1502)
                G[i+1] = 4
            elseif(N[i]<=960)
                G[i+1] = 2
            else
                G[i+1] = 3
            end
            
        elseif(G[i]==4)
            if(N[i]>=1503)
                G[i+1] = 5;
            elseif(N[i]<=970)
                G[i+1] = 3
            else
                G[i+1] = 4
            end
            
        elseif(G[i]==5)
            if(N[i]>=1504)
                G[i+1] = 6;
            elseif(N[i]<=980)
                G[i+1] = 4
            else
                G[i+1] = 5
            end
            
        elseif(G[i]==6)
            if(N[i]>=1505)
                G[i+1] = 7
            elseif(N[i]<=990)
                G[i+1] = 5
            else
                G[i+1] = 6
            end
            
        elseif(G[i]==7)
            if(N[i]>=1497)
                G[i+1] = 8
            elseif(N[i]<=1000)
                G[i+1] = 6
            else
                G[i+1] = 7
            end
            
        elseif(G[i]==8)
            if(N[i]>=1489)
                G[i+1] = 9
            elseif(N[i]<=1006)
                G[i+1] = 7
            else
                G[i+1] = 8
            end
            
        elseif(G[i]==9)
            if(N[i]>=1481)
                G[i+1] = 10
            elseif(N[i]<=1012)
                G[i+1] = 8
            else
                G[i+1] = 9
            end
            
        elseif(G[i]==10)
            if(N[i]>=1473)
                G[i+1] = 11
            elseif(N[i]<=1018)
                G[i+1] = 9
            else
                G[i+1] = 10
            end
            
        elseif(G[i]==11)
            if(N[i]>=1465)
                G[i+1] = 12
            elseif(N[i]<=1024)
                G[i+1] = 10
            else
                G[i+1] = 11
            end 
            
        elseif(G[i]==12)
            if(N[i]<=1030)
                G[i+1] = 11
            else
                G[i+1] = 12
            end
        end 
        
        i_t[i] = i_g[G[i]]
        
        N[i+1] = (30*i_g[G[i]]*i_f*(v[i]))/(3.6*pi*r_w)
        
        c_1 = r_w/(J_w+(mass*r_w*r_w)+(eta_f*i_f*i_f*eta_t*i_t[i]*i_t[i]*J_e))
        
        c_3 = eta_t*i_t[i]*eta_f*i_f
        
        del_max[i] = (a_del*N[i]*N[i])+(b_del*N[i])+c_del
        
        e[i] = -v_ref + v[i]
        
        #Step function:
        if(e[i]>0)
            u[i] = 1
        end
        
        #Computing optimal cost to go from i to each j  
        for j in j2
            p[k,i,j] = ((v[j]/3.6)-(v[i]/3.6)-(c_1*c_3*a_e*(T)*N[i])-(c_1*c_3*c_e*(T))+(c_1*c_2*(T)*((v[i]*v[i])/(3.6*3.6)))+(c_1*c_8*(T)*sin(alpha[k]+c_7)))/(c_1*c_3*b_e*(T)*del_max[i]);
            b[k,i,j] = (-(v[j]/3.6)+(v[i]/3.6)+((T)*c_1*c_3*a_d*N[i])+((T)*c_1*c_3*b_d)-((T)*c_1*c_2*((v[i]*v[i])/(3.6*3.6)))-((T)*c_1*c_8*sin(alpha[k]+c_7)))/((T)*c_1*k_b);
            
            #Checking for inefeasibility 
            if((p[k,i,j]>1.1 || p[k,i,j]<0)||((p[k,i,j]>1.1 || p[k,i,j]<0)&&(abs(v[j]-v_max)<del)))
                if(b[k,i,j]<0||b[k,i,j]>1)
                    J[k,i,j] = 10000
                    b[k,i,j] = 0
                    p[k,i,j] = 0
                else
                    p[k,i,j] = 0
                    gr[k,i,j] = G[i]
                    mf_dot = 0
                    m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                    #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j])
                    J[k,i,j] =(q1*.1*m_f[k,i,j])+(q2*2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j])
                end
            else
                b[k,i,j] = 0
                gr[k,i,j] = G[i]
                mf_dot = c_5*N[i]*p[k,i,j]*del_max[i]
                m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j]) 
                J[k,i,j] =(q1*.1*m_f[k,i,j]) + (q2*2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j])
            end
            j = j+1
        end
        i = i+1
    end
    k = k+1
end

k_1 = collect(Int,(n-1):-1:1)
i_1 = collect(Int,1:1:101)
j_1 = collect(Int,1:1:101)
n1 = length(k_1)
m1 = length(i_1)
Ter = zeros(n1+1,m1)
#T[20,:] = 10000
#T[20,51] = 0
PH = zeros(n1,m1,m1)
O = zeros(n1,m1)
Path = zeros(n1,m1)
P = zeros(n1,m1)
B = zeros(n1,m1)
gear = zeros(n1,m1)
mf = zeros(n1,m1)
lm = 0
k1 = (n-1)
i1 = 1
j1 = 1

#With terminal cost = 0, computing the optimal path to v[1] = 80 in stage 1
for k1 in k_1
    
    for i1 in i_1
        
        for j1 in j_1
            PH[k1,i1,j1] = Ter[(k1+1),j1] + J[k1,i1,j1]
            j1 = j1+1
        end 
        
        O[k1,i1] = minimum(PH[k1,i1,:])#min cost from i
        Path[k1,i1] = indmin(PH[k1,i1,:])#transition from i to j which yields min-cost
        lm = indmin(PH[k1,i1,:])
        P[k1,i1] = p[k1,i1,lm]
        B[k1,i1] = b[k1,i1,lm]
        gear[k1,i1] = gr[k1,i1,lm]
        mf[k1,i1] = m_f[k1,i1,lm]
        Ter[k1,i1] = O[k1,i1]
        i1 = i1+1
    end
    k1 = k1-1
end

z = 2
opt = zeros((n),1)
opt[1] = Path[1,51]#Starting at v[1] = 80 kmph
opt_p = zeros((n-1),1)
opt_b = zeros((n-1),1)
opt_g = zeros((n-1),1)
opt_p[1] = P[1,opt[1]]
opt_b[1] = B[1,opt[1]]
opt_mf = zeros((n-1),1)
vel = zeros((n-1),1)
vel[1] = v[opt[1]]

for z in 1:1:(n-1)
    opt[z+1] = Path[z,opt[z]]
    opt_p[z] = P[z,opt[z]]
    opt_b[z] = B[z,opt[z]]
    opt_g[z] = gear[z,opt[z]]
    opt_mf[z] = m_f[z,opt[z]]
    vel[z] = v[opt[z]]
    z = z+1
end
fuel = sum(opt_mf[1:(n-1)])
figure(1)
xlin1 = collect(1:50:2500)
xlin = collect(1:50:2500)
plin = collect(opt_p[1:1:(n-1)])
blin = collect(opt_b[1:1:(n-1)])
xlabel("Distance(m)")
ylabel("Accelerator/Brake Input")
title("Control Input(Accelerator/Brake)")
plot(xlin1,plin,"-",label = "Accelerator")
plot(xlin1,blin,"o", label = "Brake")
legend(loc = "upper right")
axis("tight")
grid("on")

figure(2)
vlin = collect(vel[1:1:(n-1)])
ref = collect(1:1:(n-1))
ref[1:(n-1)] = 85
plot(xlin1,vlin,"-",label="MPC with Fixed Weights")
#plot(xlin,ref,"-",label="Reference Speed")
grid("on")
xlabel("Distance(m)")
ylabel("Velocity(kmph)")
title("Controller Performance with 3% Road Grade")
legend(loc = "upper right")
axis("tight")

figure(3)
alin = collect((alpha[1:1:(n-1)]))
plot(xlin,alin,"-",label="Road Grade")
grid("on")
xlabel("Distance(m)")
ylabel("Altitude(m)")
legend(loc = "upper right")
title("Road Profile")
axis("tight")

figure(4)
glin = collect((opt_g[1:1:(n-1)]))
plot(xlin,glin,"-",label="Gear Number")
grid("on")
xlabel("Distance(m)")
ylabel("Gear Position")
legend(loc = "upper right")
title("Gear Position")
axis("tight")


mean(vel)

opt_mf

sum(opt_mf)

ind = 1
fueluse = zeros((n-1),1)
for ind in 1:1:(n-2)
    fueluse[ind+1] = fueluse[ind]+opt_mf[ind+1]
end

(833-909)/909

opt_mf

sum(opt_mf)


