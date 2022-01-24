#!/usr/bin/env/python
#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#Aquí se coloca la función con las ecuaciones diferenciales, donde la salida son las derivadas
def ecu_dif(t,var,Par):
    omega=Par[0];
    Eo=Par[1];
    l=Par[2];
    q=Par[3];
    m=Par[4];
    y=var[0];
    Vy=var[1];
    z=var[2];
    Vz=var[3];
    if -l<=z and z<=l:
        E=Eo*np.cos(omega*t);
    else:
        E=0;
    #fin if 
    f0=Vy;
    f1=omega*Vz;
    f2=Vz;
    f3=-omega*Vy+(q/m)*E;
    resp=np.array([f0,f1,f2,f3]);
    return resp
#fin ecuaciones diferenciales

#Esta función halla las posiciones y velocidades después de un dt
def hallar_siguienteRK4(t,var,Par,paso):
    h=paso;
    k1=ecu_dif(t,var,Par);
    k2=ecu_dif((t+0.5*h),(var+0.5*k1*h),Par);
    k3=ecu_dif((t+0.5*h),(var+0.5*k2*h),Par);
    k4=ecu_dif((t+h),(var+k3*h),Par);
    var_siguiente=var+(1/6)*h*(k1+2*k2+2*k3+k4);
    return var_siguiente
#fin hallar siguiente con Runge Kutta 4C

#Esta función halla las posiciones y velocidades en cada tiempo menor a t
def hallar_perfil(t_fin,var0,Par,Finura):
    paso=t_fin/Finura;
    ts=np.linspace(0,t_fin,Finura);
    Var=[];
    v_act=var0;
    for i in range(Finura):
        Var.append(v_act);
        var_sig=hallar_siguienteRK4(ts[i],v_act,Par,paso);
        v_act=var_sig;
    #fin for
    Var=np.array(Var);
    return [ts,Var];
#fin hallar perfil

# ---------------------------------------------------------------------------------------------------
qp=1.6022e-19; #Coulomb  Esta es la carga del protón
mp=1.6726e-27; #Kg    Esta es la masa del protón
E_final_MeV=200; #MeV Energía final
E_final_J=1.60218e-13*E_final_MeV; #J Energía final en Joules
R=1; #m    Radio final en metros.
omega=(1/R)*(((2*E_final_J)/mp)**0.5);
# omega=2*np.pi;
n_vueltas=48;
t_fin=n_vueltas*np.pi/omega;
Datos=hallar_perfil(t_fin,np.array([0,0,0,0]),[omega,5*830000/(2*0.005),0.005,qp,mp],10001);
t=Datos[0];
VarS=Datos[1];
y=VarS[:,0];
z=VarS[:,2];
Vy=VarS[:,1];
Vz=VarS[:,3];
r=[];
V=[];
r_calc=[];
r_xy=[];
t_sp=[];
tp=[];
t_act=[t[0]];
at=0;
for i in range(len(t)):
    r_dato=((y[i]**2)+(z[i]**2))**0.5;
    v_dato=((Vy[i]**2)+(Vz[i]**2))**0.5;
    if i==0:
        pass;
    else:
        if y[i-1]*y[i]<0:
            r_x_y=r_dato;
            r_cl=v_dato/omega;
            t_s_p=t[i];
            if at==0:
                t_p=t[i];
                t_act=[t[i]];
            else:
                t_p=t[i]-t_act[at-1];
            #fin if
            t_act.append(t[i]);
            r_xy.append(r_x_y);
            r_calc.append(r_cl);
            t_sp.append(t_s_p);
            tp.append(t_p);
            at=at+1;
        #fin if
    #fin if 
    r.append(r_dato);
    V.append(v_dato);
#fin foe 
r=np.array(r);
V=np.array(V);
r_xy=np.array(r_xy);
r_calc=np.array(r_calc);
t_sp=np.array(t_sp);
tp=np.array(tp);
print("************************************************************************************");
print("El último valor de r es: ",r[-1]," m");
fig1=plt.figure(figsize=(20,20));
ax1=fig1.add_subplot(211);
ax1.plot(t,VarS[:,0],label="Posición y");
ax1.plot(t,VarS[:,2],label="Posición z");
ax2=fig1.add_subplot(212);
ax2.plot(t,VarS[:,1],label="Velocidad en y");
ax2.plot(t,VarS[:,3],label="Velocidad en z");
ax1.set_title("Posiciones en función del tiempo", fontsize='x-large');
ax1.set_xlabel("tiempo [s]", fontsize='x-large');
ax1.set_ylabel("posición [m]", fontsize='x-large');
ax1.legend(loc="upper left",shadow=True);
ax2.set_title("Velocidades en función del tiempo", fontsize='x-large');
ax2.set_xlabel("tiempo [s]", fontsize='x-large');
ax2.set_ylabel("velocidad [m/s]", fontsize='x-large');
ax2.legend(loc="upper right",shadow=True);
plt.savefig("Grácifas_posicion_tiempo.png");
fig2=plt.figure(figsize=(20,20));
ax1=fig2.add_subplot(111);
ax1.plot(VarS[:,0],VarS[:,2]);
ax1.set_title("Trayectoria de la partícula", fontsize='x-large');
ax1.set_xlabel("posición en y [m]", fontsize='x-large');
ax1.set_ylabel("posición en z [m]", fontsize='x-large');
plt.savefig("Grácifa_trayectoria.png");
fig3=plt.figure(figsize=(20,20));
ax1=fig3.add_subplot(211);
ax1.plot(y,Vy);
ax1.set_title("Posición en y en función de la velocidad en y", fontsize='x-large');
ax1.set_ylabel("velocidad en y [m/s]", fontsize='x-large');
ax1.set_xlabel("posición en y [m]", fontsize='x-large');
ax2=fig3.add_subplot(212);
ax2.plot(z,Vz);
ax2.set_title("Posición en z en función de la velocidad en z", fontsize='x-large');
ax2.set_ylabel("velocidad en z [m/s]", fontsize='x-large');
ax2.set_xlabel("posición en z [m]", fontsize='x-large');
plt.savefig("Grácifa_y_Vy.png");
fig4=plt.figure(figsize=(20,20));
ax1=fig4.add_subplot(211);
ax1.scatter(r_xy,r_calc);
ax1.set_title("r obtenido de x,y vs r obtenido de p", fontsize='x-large');
ax1.set_xlabel("r obtenido de x,y [m]", fontsize='x-large');
ax1.set_ylabel("r obtenido de p [m]", fontsize='x-large');
ax2=fig4.add_subplot(212);
ax2.scatter(t_sp,tp);
ax2.set_title("Periodo en función del tiempo", fontsize='x-large');
ax2.set_xlabel("tiempo [s]", fontsize='x-large');
ax2.set_ylabel("periodo [s]", fontsize='x-large');
plt.savefig("Graficas_precisión.png");
