import os

text="""#!/bin/bash 
.\test.exe {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} 
"""
PE = [1.0]   #Atividade
G = [5.0]    #Intensidade do acoplamento
PHI = [0.6]  #Fração de empacotamento (densidade)
N = [200]    #Número de aneis
Ka = [1.0]   #Intensidade da força de área
Kl = [40.0] #Intensidade da Constante de mola do anel
Kc = [40.0] #Intensidade da força de repulsão
Kadh = [1.0] #Intensidade da força de adesão
p0 = [3.8]   #Tensão cortical
radh = [1.5] #Alcance da força
par=[]

for i1 in PE:
    for i2 in G:
        for i3 in PHI:
            for i4 in N:
                for i5 in Ka:
                    for i6 in Kl:         
                        for i7 in Kc:
                            for i8 in Kadh:
                                for i9 in p0:
                                    for i10 in radh:
                                        par.append([i1,i2,i3,i4,i5,i6,i7,i8,i9,i10])	 																			
                              
for i in par:
	
    o=open('torun_nocluster.sh','w')
    o.write(text.format(*i))
    #o.write("#!/bin/bash \n")
    #o.write("#SBATCH -n 1 # Number of cores\n")
    #o.write("#SBATCH -N 1 # Number of nodes\n")
    #o.write("#SBATCH -t 03-00:00 # Tempo limite de execucao (D-HH:MM)\n")
    #o.write("#SBATCH -p short # Partition to submit to \n")
    #o.write("#SBATCH --qos qos_short # QOS\n")
    #o.write("./a.out %f %f %f %f %f %f %f\n"%(i[0],i[1],i[2],i[3],i[4],i[5],i[7]))
    o.close()	
    #os.system('sbatch torun.scpt')
    os.system('sh torun_nocluster.sh')





