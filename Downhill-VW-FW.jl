
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
m_f[1,:,:] = 100
J = zeros(n,m,m)
del_max = zeros(m,1)#Maximum Fueling
e = zeros(m,1)
i_t = zeros(m,1)
N = zeros(m+1,1)
G = zeros(m+1,1)
#Road slope
alpha = zeros(n,1)
alpha[15:25] = -0.03
u = zeros(m,1)

#Intial Conditions
N[1] = (30*minimum(v)*i_f)/(3.6*pi*r_w)
G[1] = 12

#Cost Function Parameters:
q1=zeros(n,1)
q2=zeros(n,1)
q3=zeros(n,1)
q4=zeros(n,1)
q5=zeros(n,1)
del = 0.1
i = 1
j = 1
k = 1

#Stage
for k in k2
    if(alpha[k]==0)
            q1[k] = 1
            q2[k] = 5
            q3[k] = 15
            q5[k] = 0
    elseif(alpha[k]>0)
            q1[k] = 0.1
            q2[k] = 1
            q3[k] = 15
            q5[k] = 0
    else
            q1[k] = 25
            q2[k] = 0.001
            q3[k] = 15
            q5[k] = 100
    end
    
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
                G[i+1] = 5
            elseif(N[i]<=970)
                G[i+1] = 3
            else
                G[i+1] = 4
            end
            
        elseif(G[i]==5)
            if(N[i]>=1504)
                G[i+1] = 6
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
        
        e[i] = abs(-v_ref + v[i])
        
        #Step function:
        if(e[i]>0)
            u[i] = 1
        end

        
        #Computing optimal cost to go from i to each j  
        for j in j2
            p[k,i,j] = ((v[j]/3.6)-(v[i]/3.6)-(c_1*c_3*a_e*(T)*N[i])-(c_1*c_3*c_e*(T))+(c_1*c_2*(T)*((v[i]*v[i])/(3.6*3.6)))+(c_1*c_8*(T)*sin(alpha[k]+c_7)))/(c_1*c_3*b_e*(T)*del_max[i]);
            b[k,i,j] = (-(v[j]/3.6)+(v[i]/3.6)+((T)*c_1*c_3*a_d*N[i])+((T)*c_1*c_3*b_d)-((T)*c_1*c_2*((v[i]*v[i])/(3.6*3.6)))-((T)*c_1*c_8*sin(alpha[k]+c_7)))/((T)*c_1*k_b);
            
            #Checking for inefeasibility 
            if((p[k,i,j]>1.1 || p[k,i,j]<0.1)||((p[k,i,j]>1.1 || p[k,i,j]<0.1)&&(abs(v[j]-v_max)<del)))
                if(b[k,i,j]>0.1||b[k,i,j]<-1.1)
                    J[k,i,j] = 10000
                    p[k,i,j] = 0
                    b[k,i,j] = 0
                else
                    p[k,i,j] = 0
                    if(b[k,i,j]>0)
                        b[k,i,j] = 0
                    end
                    gr[k,i,j] = G[i]
                    mf_dot = 0
                    m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                    #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j])
                    J[k,i,j] =(q1[k]*.01*m_f[k,i,j])+(q2[k]*2*u[i]*e[i]*e[i])+(q3[k]*abs(-v[j]+v[i]))+(q5[k]*abs(b[k,i,j]))
                end
            else
                b[k,i,j] = 0
                gr[k,i,j] = G[i]
                mf_dot = c_5*N[i]*(p[k,i,j])*del_max[i]
                m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j]) 
                J[k,i,j] =(q1[k]*.01*m_f[k,i,j]) + (q2[k]*2*u[i]*e[i]*e[i])+(q3[k]*abs(-v[j]+v[i]))+(q5[k]*abs(b[k,i,j]))
            end      
        end
    end
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
mf = zeros(n1,m1)
lm = 0
k1 = (n-1)
i1 = 1
j1 = 1
gear = zeros(n1,m1)

#With terminal cost = 0, computing the optimal path to v[1] = 80 in stage 1
for k1 in k_1
    
    for i1 in i_1
        
        for j1 in j_1
            PH[k1,i1,j1] = Ter[(k1+1),j1] + J[k1,i1,j1]
        end 
        
        O[k1,i1] = minimum(PH[k1,i1,:])#min cost from i
        Path[k1,i1] = indmin(PH[k1,i1,:])#transition from i to j which yields min-cost
        lm = indmin(PH[k1,i1,:])
        P[k1,i1] = p[k1,i1,lm]
        B[k1,i1] = b[k1,i1,lm]
        mf[k1,i1] = m_f[k1,i1,lm]
        gear[k1,i1] = gr[k1,i1,lm]
        Ter[k1,i1] = O[k1,i1]
    end
end
fuel = 0
z = 1
opt = zeros((n),1)
opt[1] = Path[1,51]#Starting at v[1] = 80 kmph
opt_p = zeros((n-1),1)
opt_b = zeros((n-1),1)
opt_mf = zeros((n-1),1)
opt_p[1] = P[1,opt[1]]
opt_b[1] = B[1,opt[1]]
opt_g = zeros((n-1),1)
vel = zeros((n-1),1)
vel[1] = v[opt[1]]

for z in 1:1:(n-1)
    opt[z+1] = Path[z,opt[z]]
    opt_p[z] = P[z,opt[z]]
    opt_b[z] = B[z,opt[z]]
    opt_mf[z] = mf[z,opt[z]]
    opt_g[z] = gear[z,opt[z]]
    vel[z] = v[opt[z]]
end
fuel=sum(opt_mf[1:(n-1)])

figure(1)
xlin = collect(1:50:2500)
plin = collect(opt_p[1:1:(n-1)])
blin = collect(opt_b[1:1:(n-1)])
xlabel("Distance(m)")
ylabel("Accelerator/Brake Input")
title("Control Input(Accelerator/Brake)")
plot(xlin,plin,"-",label = "Accelerator")
plot(xlin,blin,"-", label = "Brake")
legend(loc = "lower right")
axis("tight")
grid("on")

figure(2)
vlin = collect(vel[1:1:(n-1)])
ref = collect(1:1:(n-1))
ref[1:(n-1)] = 85
plot(xlin,vlin,"-",label="MPC")
#plot(xlin,ref,"-",label="Reference Speed")
grid("on")
xlabel("Distance(m)")
ylabel("Velocity(kmph)")
title("Controller Performance with -3% Road Grade")
legend(loc = "lower right")
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



fuel

opt_mf

maximum(opt_mf)

fuel

fig = figure("Performance Of Optimal & PI Control",figsize=(15,15))
subplots_adjust(hspace=0.3) # Set the vertical spacing between axes

subplot(411) # Create the 1st axis of a 3x1 array of axes
ax1 = gca()
alt = zeros((n-1),1)
alt[1:15] = 0.03*xlin[15]
malt=1
for malt in 1:1:10
    alt[(15+malt)] = 0.03*xlin[(15+malt)]
end
alt[25:(n-1)] = 0.03*xlin[25] 

xticks(0:500:2500,fontsize = 20) # Set the x axis to a logarithmic scale
xlim(0,2500)
yticks(42:-7:21,fontsize = 20)
plot(xlin,alt, color = "b","-",color = "brown",linewidth = 4,label = "Road Profile")# Set the y-tick range and step size, 0.1 to 0.9 in increments of 0.2
ylim(42,20) # Set the y-limits from 0.0 to 1.0
xlabel("Distance(m)",fontsize = 20)
ylabel("Altitude(m)",fontsize = 20)
title("Road Grade = -3%",fontsize = 20)
legend(loc="upper right",fontsize = 20)
grid("on")
#axis("tight")

subplot(412) # Create the 1st axis of a 3x1 array of axes
ax2 = gca()
xticks(0:500:2500,fontsize = 20) # Set the x axis to a logarithmic scale
xlim(0,2500)
yticks(80:5:90,fontsize = 20) # Set the y-tick range and step size, 0.1 to 0.9 in increments of 0.2
ylim(80,90)
plot(xlin,vlin,"-",label="Velocity(Varying-Weight)",linewidth = 4)
plot(xlin,ref,"_",color="r",label="Reference Speed",linewidth = 4)
#plot(xlin,vlin1,"_",color="b",label = "PI")
plot(xlin,vlin2,"--",color ="black",label = "Velocity(Fixed-Weight)",linewidth = 4)
annotate("Decline in speed",fontsize = 20,
xy=[2250;81],# Arrow tip
    xytext=[1500;82], # Text offset from tip
    xycoords="data", # Coordinates in in "data" units
    arrowprops=["facecolor"=>"black"])
xlabel("Distance(m)",fontsize = 20)
ylabel("Velocity(kmh)",fontsize = 20)
legend(loc = "upper right",fontsize = 15)
grid("on")
#axis("tight")

subplot(413) # Create the 2nd axis of a 3x1 array of axes
ax3 = gca()
#setp(ax2[:get_xticklabels](),visible=false)
xticks(0:500:2500,fontsize = 20) # Set the x axis to a logarithmic scale
xlim(0,2500)
yticks(0:0.5:1,fontsize = 20)
ylim(0.0,1.1)
plot(xlin,plin,"-",label="Accelerator(Varying-Weight)",linewidth = 4)
plot(xlin,blin,"o",label="Brake(Varying-Weight)",linewidth = 4)
plot(xlin,plin2,"--",color = "black",label="Accelerator(Fixed-Weight)",linewidth = 4)
plot(xlin,blin2,"-",color = "black",label="Brake(Fixed-Weight)",linewidth = 4)
xlabel("Distance(m)",fontsize = 20)
ylabel("Control Inputs",fontsize = 20)
legend(loc = "upper right",fontsize = 15)
grid("on")
#axis("tight")

subplot(414) # Create the 3rd axis of a 3x1 array of axes
xticks(0:500:2500,fontsize = 20) # Set the x axis to a logarithmic scale
xlim(0,2500)
xlabel("Distance(m)")
grid("on")
yticks(1:1:12,fontsize = 15)
ylim(1,13)
plot(xlin,glin,"-",label="Gear(Varying-Weight)",linewidth = 4)
#plot(xlin,glin1,"_",label="PI")
plot(xlin,glin2,"_",color = "black",label="Gear(Fixed-Weight)",linewidth = 4)
ylabel("Gear Position",fontsize = 20)
xlabel("Distance(m)",fontsize = 20)
legend(loc ="lower right",fontsize = 20)
#axis("tight")
fig[:canvas][:draw]() # Update the figure
suptitle("Performance Analyses Of Fixed-Weight Optimal Control Vs Varying-Weight Optimal Control",fontsize = 20)

opt_mf

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
m_f[1,:,:] = 100
J = zeros(n,m,m)
del_max = zeros(m,1)#Maximum Fueling
e = zeros(m,1)
i_t = zeros(m,1)
N = zeros(m+1,1)
G = zeros(m+1,1)
#Road slope
alpha = zeros(n,1)
alpha[15:25] = -0.03
u = zeros(m,1)

#Intial Conditions
N[1] = (30*minimum(v)*i_f)/(3.6*pi*r_w)
G[1] = 12

#Cost Function Parameters:
q1=2
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
            if((p[k,i,j]>1.1 || p[k,i,j]<-0.1)||((p[k,i,j]>1.1 || p[k,i,j]<-0.1)&&(abs(v[j]-v_max)<del)))
                if(b[k,i,j]>0.1||b[k,i,j]<-1.1)
                    J[k,i,j] = 10000
                    b[k,i,j] = 0
                    p[k,i,j] = 0
                else
                    p[k,i,j] = 0
                    gr[k,i,j] = G[i]
                    mf_dot = 0
                    m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                    #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j])
                    J[k,i,j] =(q1*0.1*m_f[k,i,j])+(q2*2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*abs(b[k,i,j]))
                end
            else
                b[k,i,j] = 0
                gr[k,i,j] = G[i]
                mf_dot = c_5*N[i]*abs(p[k,i,j])*del_max[i]
                m_f[k+1,i,j] = m_f[k,i,j]+ (T)*mf_dot;
                #J[k,i,j] = (q1*m_f[k,i,j])+(q2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*b[k,i,j]) 
                J[k,i,j] =(q1*0.1*m_f[k,i,j]) + (q2*2*u[i]*e[i]*e[i])+(q3*abs(-v[j]+v[i]))+(q5*abs(b[k,i,j]))
            end      
        end
    end
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
mf = zeros(n1,m1)
lm = 0
k1 = (n-1)
i1 = 1
j1 = 1
gear = zeros(n1,m1)

#With terminal cost = 0, computing the optimal path to v[1] = 80 in stage 1
for k1 in k_1
    
    for i1 in i_1
        
        for j1 in j_1
            PH[k1,i1,j1] = Ter[(k1+1),j1] + J[k1,i1,j1]
        end 
        
        O[k1,i1] = minimum(PH[k1,i1,:])#min cost from i
        Path[k1,i1] = indmin(PH[k1,i1,:])#transition from i to j which yields min-cost
        lm = indmin(PH[k1,i1,:])
        P[k1,i1] = p[k1,i1,lm]
        B[k1,i1] = b[k1,i1,lm]
        mf[k1,i1] = m_f[k1,i1,lm]
        gear[k1,i1] = gr[k1,i1,lm]
        Ter[k1,i1] = O[k1,i1]
        if(P[k1,i1]>=1)
            P[k1,i1] = 1
        elseif(P[k1,i1]<=0)
            P[k1,i1] = 0
        end
        if(B[k1,i1]>=0)
            B[k1,i1] = 0
        elseif(B[k1,i1]<=-1)
            B[k1,i1] = -1
        end
    end
end

z = 2
opt = zeros((n),1)
opt[1] = Path[1,51]#Starting at v[1] = 80 kmph
opt_p = zeros((n-1),1)
opt_b = zeros((n-1),1)
opt_mf = zeros((n-1),1)
opt_p[1] = P[1,opt[1]]
opt_b[1] = B[1,opt[1]]
opt_g = zeros((n-1),1)
vel = zeros((n-1),1)
vel[1] = v[opt[1]]

for z in 1:1:(n-1)
    opt[z+1] = Path[z,opt[z]]
    opt_p[z] = P[z,opt[z]]
    opt_b[z] = B[z,opt[z]]
    opt_mf[z] = mf[z,opt[z]]
    opt_g[z] = gear[z,opt[z]]
    vel[z] = v[opt[z]]
end
plin2 = collect(opt_p[1:1:(n-1)])
blin2 = collect(opt_b[1:1:(n-1)])
vlin2 = collect(vel[1:1:(n-1)])
glin2 = collect((opt_g[1:1:(n-1)]))

fuel

alt

1751-201

xlin[25]

1201*0.03



opt_mf


