import math
import csv
import matplotlib.pyplot as plt

def NINT(x):
    return int(math.copysign(math.floor(abs(x) + 0.5), x))


DR = 0.01
D_O = 3.05
D_S = 4.35
D_T = 3.62
FMAX = 0.01
FP = 1.0
ICODE = 0
NP = 805
PI = 3.141592653589793
RBT = 8.0
TOLE = 0.0001

RHOBS = 0.0304
RHOBO = 2.2619
RHOBT = 2.2923

RHOB = [RHOBS * 0.00060221415,RHOBO * 0.00060221415,RHOBT * 0.00060221415]
DP = [D_S, D_O, D_T]

IR = [0] * 3
NR = [0] * 3

for J in range (3):
    IR[J]=NINT(DP[J]/2./DR)+1
    NR[J]=NINT((RBT-DP[J]/2.)/DR)+1

IRMIN = min(IR)
NRMAX = max(NR)

DD = [[[0.0 for _ in range(NP + 1)] for _ in range(6)] for _ in range(3)]

def WD(I,IR,DP,RHO,NI,NR):

    R = float(I-1) * DR

    n0i = [0.0] * 3
    n1i = [0.0] * 3
    n2i = [0.0] * 3
    n3i = [0.0] * 3
    nv1i = [0.0] * 3
    nv2i = [0.0] * 3



    for J in range(0, 3):
        ...
        RIN=R-DP[J]/2.
        RIP=R+DP[J]/2.
        KIN=NINT(RIN/DR)+1
        KIP=NINT(RIP/DR)+1

        if KIN<=IR[J]:
            KIN=IR[J]

        if KIP>= NINT((RBT-DP[J]/2.)/DR)+1: 
            KIP=NINT((RBT-DP[J]/2.)/DR)+1 

        if KIN<KIP:
            n2i[J]=0.0
            n3i[J]=0.0
            nv2i[J]=0.0
            # TODO: review loop bounds
            for K in range(KIN, KIP+1):
                R1=float(K-1)*DR
                if (I==KIN) or (I==KIP):
                    n2i[J]=n2i[J]+RHO[J][K]/2.
                    n3i[J]=n3i[J]+RHO[J][K]*(DP[J]**2/4.-(R1-R)**2)/2.
                    nv2i[J]=nv2i[J]+RHO[J][K]*(R1-R)/2.
                else:
                    n2i[J]=n2i[J]+RHO[J][K]
                    n3i[J]=n3i[J]+RHO[J][K]*(DP[J]**2/4.-(R1-R)**2)
                    nv2i[J]=nv2i[J]+RHO[J][K]*(R1-R)
            
            n2i[J]=n2i[J]*PI*DP[J]*DR
            n3i[J]=n3i[J]*PI*DR
            nv2i[J]=nv2i[J]*2.*PI*DR
            n0i[J]=n2i[J]/(PI*DP[J]**2)
            n1i[J]=n2i[J]/(2.*PI*DP[J])
            nv1i[J]=nv2i[J]/(2.*PI*DP[J])
            NI[J][0]=n0i[J]
            NI[J][1]=n1i[J]
            NI[J][2]=n2i[J]
            NI[J][3]=n3i[J]
            NI[J][4]=nv1i[J]
            NI[J][5]=nv2i[J]
        else:
            NI[J][0]=0.0
            NI[J][1]=0.0
            NI[J][2]=0.0
            NI[J][3]=0.0
            NI[J][4]=0.0
            NI[J][5]=0.0



def EXTERN(K,IR,NR,EXT):

    for J in range(0, 3):
        if K<IR[J] or K>NR[J]:
            EXT[J]=1.e50
        else:
            EXT[J]=0.0


def DFEXRC(K,DP,DD,DFRC,IR,NR):

    R=float(K-1)*DR

    for J in range(3):
        RIN=R-DP[J]/2.0
        RIP=R+DP[J]/2.0
        KIN=NINT(RIN/DR)+1
        KIP=NINT(RIP/DR)+1

        if KIN<=IR[J]:
            KIN=IR[J]

        if KIP>= NINT((RBT-DP[J]/2.)/DR)+1: 
            KIP=NINT((RBT-DP[J]/2.)/DR)+1 

        DFRC1 = [0.0] * 6

        # TODO: review loop bounds
        for I in range(KIN, KIP+1):
            R1 = float(I-1) * DR
            if (I==KIN) or (I==KIP):
                DFRC1[0]=DFRC1[0]+DD[J][0][I-1]/2.
                DFRC1[1]=DFRC1[1]+DD[J][1][I-1]/2.
                DFRC1[2]=DFRC1[2]+DD[J][2][I-1]/2.
                DFRC1[3]=DFRC1[3]+DD[J][3][I-1]*(DP[J]**2/4.-(R1-R)**2)/2.
                DFRC1[4]=DFRC1[4]+DD[J][4][I-1]*(R-R1)/2.
                DFRC1[5]=DFRC1[5]+DD[J][5][I-1]*(R-R1)/2.
            else:
                DFRC1[0]=DFRC1[0]+DD[J][0][I-1]
                DFRC1[1]=DFRC1[1]+DD[J][1][I-1]
                DFRC1[2]=DFRC1[2]+DD[J][2][I-1]
                DFRC1[3]=DFRC1[3]+DD[J][3][I-1]*(DP[J]**2/4.-(R1-R)**2)
                DFRC1[4]=DFRC1[4]+DD[J][4][I-1]*(R-R1)
                DFRC1[5]=DFRC1[5]+DD[J][5][I-1]*(R-R1)


        DFRC1[0]=DFRC1[0]*DR/DP[J]
        DFRC1[1]=DFRC1[1]*DR/2.0
        DFRC1[2]=DFRC1[2]*DR*PI*DP[J]
        DFRC1[3]=DFRC1[3]*DR*PI
        DFRC1[4]=DFRC1[4]*DR/DP[J]
        DFRC1[5]=DFRC1[5]*DR*PI*2.0

        
        DFRC[J]=0.0
        for J1 in range(0, 6):
            DFRC[J]=DFRC[J]+DFRC1[J1]


def CP(RHOB,DP,MU):

    an0i = [0.0]*3
    an1i = [0.0]*3
    an2i = [0.0]*3
    an3i = [0.0]*3

    MUID  = [0.0 for _ in range(3)]
    MUHS  = [0.0 for _ in range(3)]


    for J in range(0, 3):
        an0i[J]=RHOB[J]
        an1i[J]=RHOB[J]*DP[J]/2.0
        an2i[J]=RHOB[J]*DP[J]**2*PI
        an3i[J]=RHOB[J]*DP[J]**3*PI/6.0
 
    AN0=an0i[0]+an0i[1]+an0i[2]
    AN1=an1i[0]+an1i[1]+an1i[2]
    AN2=an2i[0]+an2i[1]+an2i[2]
    AN3=an3i[0]+an3i[1]+an3i[2]

    DHN0I=-math.log(1.-AN3)
    DHN1I=AN2/(1.-AN3)
    DHN2I=AN1/(1.-AN3)+AN2**2/(12.*PI)*(math.log(1.-AN3)/AN3**2+1./(AN3*(1.-AN3)**2))
    DHN3I=AN0/(1.-AN3)+AN1*AN2/(1.-AN3)**2-AN2**3/(36.*PI)*(2.*math.log(1.-AN3)/AN3**3+(2.-5.*AN3+AN3**2)/(AN3**2*(1.-AN3)**3))

    for I in range(3):

        if RHOB[I]>1.e-6:
            MUID[I]=math.log(RHOB[I])
        else:
            MUID[I]=-1.0e20

    # # end do

    for I in range(3):
        MUHS[I]=DHN0I+0.5*DP[I]*DHN1I+PI*DP[I]**2*DHN2I+PI*DP[I]**3/6.*DHN3I


    for I in range(3):
        MU[I]=MUID[I]+MUHS[I]



def DDRC(IR,NR,DP,RHO,DD):


    NRMAX2=max(NR[0]+NINT(0.5*DP[0]/DR), NR[1]+NINT(0.5*DP[1]/DR), NR[2]+NINT(0.5*DP[2]/DR))

    n0i = [0.0] * 3
    n1i = [0.0] * 3
    n2i = [0.0] * 3
    n3i = [0.0] * 3
    nv1i = [0.0] * 3
    nv2i = [0.0] * 3

    DDR = [0.0] * 6

    NI = [[0.0 for _ in range(6)] for _ in range(3)]

    for K in range( 1, NRMAX2 + 1):
        WD(K,IR,DP,RHO,NI,NR)
        for J in range(3):
            n0i[J]=NI[J][0]
            n1i[J]=NI[J][1]
            n2i[J]=NI[J][2]
            n3i[J]=NI[J][3]
            nv1i[J]=NI[J][4]
            nv2i[J]=NI[J][5]
        # # end do

        N0=n0i[0]+n0i[1]+n0i[2]
        N1=n1i[0]+n1i[1]+n1i[2]
        N2=n2i[0]+n2i[1]+n2i[2]
        N3=n3i[0]+n3i[1]+n3i[2]
        NV1=nv1i[0]+nv1i[1]+nv1i[2]
        NV2=nv2i[0]+nv2i[1]+nv2i[2]


        if N3<=1.0e-8:
            for J1 in range(3):
                for J2 in range(6):
                    DD[J1][J2][K]=0.0
    
        else:

            DDR[0]=-math.log(1.-N3)
            DDR[1]=N2/(1.-N3)
            DDR[2]=N1/(1.-N3)+(math.log(1.-N3)/N3+1./(1.-N3)**2)*(N2**2-NV2**2)/(12.*PI*N3)
            DDR[3]=N0/(1.-N3)+(N1*N2-NV1*NV2)/(1-N3)**2-(math.log(1.-N3)/(18.*PI*N3**3)+1./(36.*PI*N3**2*(1.-N3))+(1.-3.*N3)/(36.*PI*N3**2*(1.-N3)**3))*(N2**3-3.*N2*NV2**2)
            DDR[4]=-NV2/(1.-N3)
            DDR[5]=-NV1/(1.-N3)-(math.log(1.-N3)/N3+1./(1.-N3)**2)*N2*NV2/(6.*PI*N3)

        for j in range(3):
            DD[j][0][K]=DDR[0]
            DD[j][1][K]=DDR[1]
            DD[j][2][K]=DDR[2]
            DD[j][3][K]=DDR[3]
            DD[j][4][K]=DDR[4]
            DD[j][5][K]=DDR[5]




def DPS(IR,NR,DP,MU,DD,RHO,RHO1,RHOB):

    # (patched) do not reassign RHO1 here
    EXT  = [0.0]*3
    DFRC = [0.0]*3
    
    

    for K in range(IRMIN, NRMAX+1):
        R=float(K-1)*DR
        EXTERN(K,IR,NR,EXT)
        DFEXRC(K,DP,DD,DFRC,IR,NR)
        if K>=IR[0]:
            RHO1[0][K]=math.exp(MU[0]-EXT[0]-DFRC[0])

        if K>=IR[1]:
            RHO1[1][K]=math.exp(MU[1]-EXT[1]-DFRC[1])

        if K>=IR[2]:
            RHO1[2][K]=math.exp(MU[2]-EXT[2]-DFRC[2])
    

def main():
    print ("Program Start")
 


    RHO = [[0.0 for _ in range(NP + 1)] for _ in range(3)]
    RHO1 = [[0.0 for _ in range(NP + 1)] for _ in range(3)]
    MU = [0.0]*3
    DD = [[[0.0 for _ in range(NP + 1)] for _ in range(6)] for _ in range(3)]


    if ICODE==0:
        for J in range(0, 3):
            for K in range(0, NP):
                if K<IR[J] or K>NR[J]:
                    RHO[J][K]=0.0
                else:
                    RHO[J][K]=RHOB[J]

    if ICODE!=0:
        f1 = open("RHO", "r")
        for K in range(0, NP):
            R, RP1, RP2, RP3 = map(float, f1.readline().split())
            RHO[0][K] = RP1
            RHO[1][K] = RP2
            RHO[2][K] = RP3
        f1.close()

    CP(RHOB,DP,MU)

    ITER=1

    while True:
        DDRC (IR,NR,DP,RHO,DD)

        if (ITER!=0) and (ITER % 100 == 0):
            f11 = open("RHO", "w")
            for K in range(0, NP+1):
                R=float(K)*DR
                print(R, *(RHO[J][K] for J in range(3)), file=f11)
            f11.close()

        DPS (IR,NR,DP,MU,DD,RHO,RHO1,RHOB)

        ACCU=0
        ERR=0.0
        ACCERR=0.0
        for J in range(3):
            if RHOB[J]>1.e-6:
                for K in range(IR[J], NR[J]+1):
                    delta = abs(RHO1[J][K] - RHO[J][K])
                    if delta >= TOLE:
                        ACCU=ACCU+1
                    if ERR<delta:
                        ERR=delta
                    ACCERR=ACCERR+delta


        if ACCU>0:
            ITER=ITER+1
            if ITER % 10 == 0:
                print(ITER,ACCU,ERR,ACCERR)

            MIX_F = min(FP / ERR, FMAX)
            
            for J in range(3):
                for K in range(IR[J], NR[J]+1):
                    RHO[J][K]=(1.-MIX_F)*RHO[J][K]+MIX_F*RHO1[J][K]
            continue
        else:
             break




    with open("Reduced Density Profile.csv", "w", newline="") as f42:
        writer = csv.writer(f42)
        writer.writerow(["R (distance)", "Species1_ReducedDensity", "Species2_ReducedDensity", "Species3_ReducedDensity"])
        for K in range(0, NP):
            R = float(K) * DR
            row = [R] + [RHO[J][K] / RHOB[J] for J in range(3)]
            writer.writerow(row)

    print("✅ Reduced Density Profile.csv created. This file can be opened directly in Excel without formatting issues.")

    # Write Concentration Profile in CSV format
    with open("Concentration Profile.csv", "w", newline="") as f41:
        writer = csv.writer(f41)
        writer.writerow(["R (distance)", "Species1_Concentration", "Species2_Concentration", "Species3_Concentration"])
        for K in range(0, NP):
            R = float(K) * DR
            row = [R] + [RHO[J][K] / 0.00060221415 for J in range(3)]
            writer.writerow(row)

    print("✅ Concentration Profile.csv created. You can open it directly in Excel — columns and headers are automatically aligned.")


    r_vals = [float(K) * DR for K in range(NP)]

    # Concentration Profile Plot
    plt.figure()
    for j in range(3):
        plt.plot(r_vals, [RHO[j][K] / 0.00060221415 for K in range(NP)], label=f"Species {j + 1}")
    plt.title("Concentration Profile")
    plt.xlabel("R (distance)")
    plt.ylabel("Concentration")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Concentration_Profile.png", dpi=300)
    plt.close()
    print("📈 Concentration_Profile.png saved.")

    # Reduced Density Profile Plot
    plt.figure()
    for j in range(3):
        plt.plot(r_vals, [RHO[j][K] / RHOB[j] for K in range(NP)], label=f"Species {j + 1}")
    plt.title("Reduced Density Profile")
    plt.xlabel("R (distance)")
    plt.ylabel("Reduced Density")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Reduced_Density_Profile.png", dpi=300)
    plt.close()
    print("📈 Reduced_Density_Profile.png saved.")


    for J in range(0, 3):
        RHOB[J]=RHOB[J]/0.00060221415

    RO = [0.0] * 3

    for J in range(0, 3):
        RO[J]=0.0
        for K in range(IRMIN, NRMAX+1):
            if (K==IR[J])or(K==NR[J]):
                RO[J]=RO[J]+RHO[J][K]*DR/2.0
            else:
                RO[J]=RO[J]+RHO[J][K]*DR
        RO[J]=RO[J]/RBT/0.00060221415
    # # end do

    f60 = open("Output.dat", "w")
    print("Cation 1 Average Concentration in Pore:", RO[0], file=f60)
    print("Cation 2 Average Concentration in Pore:", RO[1], file=f60)
    print("ANIon Average Concentration in Pore:", RO[2], file=f60)
    f60.close()

    print ("Program End")
if __name__ == "__main__":
    main()