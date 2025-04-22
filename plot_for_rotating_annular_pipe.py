import matplotlib.pyplot as plt
import numpy as np                                                            
import matplotlib as mpl            
import math


def velocity_profile():
    y_plus_Om03_Ma08 = [0]*100
    y_plus_Om0_Ma08 = [0]*100
    y_plus_Om03_Ma03 = [0]*100
    y_plus_Om0_Ma03 = [0]*100
    y_plus_Om03_Ma15 = [0]*100
    y_plus_Om0_Ma15 = [0]*100
    u_tau_Om03_Ma08 = [0]*100
    u_tau_Om03_Ma03 = [0]*100
    u_tau_Om0_Ma08 = [0]*100
    u_tau_Om0_Ma03 = [0]*100
    u_tau_Om03_Ma15 = [0]*100
    u_tau_Om0_Ma15 = [0]*100
    yfit = [0.1,100]
    utau_fit = [1/0.41*math.log(y)+5.2 for y in yfit]

    fopen = open("./velocity_profile/yplus_Om0.3_Ma0.8.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om03_Ma08[i] = float(line[0])
        u_tau_Om03_Ma08[i] = float(line[1])
        i = i + 1
    fopen.close()

    fopen = open("./velocity_profile/yplus_Om0_Ma0.8.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om0_Ma08[i] = float(line[0])
        u_tau_Om0_Ma08[i] = float(line[1])
        i = i + 1
    fopen.close()

    fopen = open("./velocity_profile/yplus_Om0.3_Ma0.3.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om03_Ma03[i] = float(line[0])
        u_tau_Om03_Ma03[i] = float(line[1])
        i = i + 1
    fopen.close()

    fopen = open("./velocity_profile/yplus_Om0_Ma0.3.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om0_Ma03[i] = float(line[0])
        u_tau_Om0_Ma03[i] = float(line[1])
        i = i + 1
    fopen.close()

    fopen = open("./velocity_profile/yplus_Om0.3_Ma1.5.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om03_Ma15[i] = float(line[0])
        u_tau_Om03_Ma15[i] = float(line[1])
        i = i + 1
    fopen.close()

    fopen = open("./velocity_profile/yplus_Om0_Ma1.5.dat",'r')
    i = 1
    while i < 100:
        line_temp = fopen.readline()
        line=line_temp.split()
        y_plus_Om0_Ma15[i] = float(line[0])
        u_tau_Om0_Ma15[i] = float(line[1])
        i = i + 1
    fopen.close()

    fig = plt.figure()
    fig.set_size_inches(12, 27/4)
    # plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    ax = fig.add_subplot(111)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    plt.tick_params(labelsize=15)
    plt.tick_params(axis='both',width=2,length=6)
    plt.plot(y_plus_Om03_Ma08,u_tau_Om03_Ma08,c='royalblue',label='$\omega=0.3,Ma=0.8$',linewidth=2.5)
    plt.plot(y_plus_Om0_Ma08,u_tau_Om0_Ma08,c='darkblue',label='$\omega=0,Ma=0.8$',linewidth=2.5)
    plt.plot(y_plus_Om03_Ma03,u_tau_Om03_Ma03,c='orange',label='$\omega=0.3,Ma=0.3$',linewidth=2.5)
    plt.plot(y_plus_Om0_Ma03,u_tau_Om0_Ma03,c='peru',label='$\omega=0,Ma=0.3$',linewidth=2.5)
    plt.plot(y_plus_Om03_Ma15,u_tau_Om03_Ma15,c='limegreen',label='$\omega=0.3,Ma=1.5$',linewidth=2.5)
    plt.plot(y_plus_Om0_Ma15,u_tau_Om0_Ma15,c='forestgreen',label='$\omega=0,Ma=1.5$',linewidth=2.5)
    plt.plot(yfit,utau_fit,'k',label='$1/0.41 ln(y^+)+5.1$',linewidth=2.5,linestyle='--')
    plt.xlabel('$y^+$', rotation=0, fontdict={'size':20}, labelpad=-0.1)
    plt.xscale('log')
    # plt.yscale('log')
    plt.ylabel('$u$', rotation=0, fontdict={'size':20}, labelpad=-0.3)
    plt.xlim(0.1,80)
    plt.ylim(0,20)
    plt.legend(loc='upper left', fontsize=20, frameon=False)
    plt.savefig('./velocity_profile/velocity_profile.png',dpi=200,pad_inches=0.1)
    plt.close()

velocity_profile()