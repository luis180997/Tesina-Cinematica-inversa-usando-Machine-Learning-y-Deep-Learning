import numpy as np
from numpy import sin, cos, pi

def matriz_homogenea_DH(theta,d,a,alpha):
    
    ct = cos(theta)
    st = sin(theta)
    ca = cos(alpha)
    sa = sin(alpha)
    
    # VERIFICAMOS EL ANGULO "theta"
    if(abs(theta)==pi or abs(theta) == 2*pi):
        st = 0
    elif(abs(theta) == pi/2):
        ct = 0
        
    # VERIFICAMOS EL ANGULO "alpha"
    if(abs(alpha)==pi or abs(alpha) == 2*pi):
        sa = 0
    elif(abs(alpha) == pi/2):
        ca = 0
    
    h = [[ct, -st*ca, st*sa, a*ct],
         [st, ct*ca, -ct*sa, a*st],
         [0, sa, ca, d],
         [0, 0, 0, 1]]
    H = np.array(h)
    
    return H

def angulos_euler(H, flag1=1):

    if(np.shape(H) != (3,3)):
        R = H[0:3,0:3]
    else:
        R = H

    if (np.abs(R[2,2] != 1)): # R[0,2] != 0 and R[1,2] != 0
        theta = np.arctan2(flag1*np.sqrt(1-R[2,2]**2),R[2,2])
        phi = np.arctan2(flag1*R[1,2],flag1*R[0,2])
        psi = np.arctan2(flag1*R[2,1],-flag1*R[2,0])

    elif(np.abs(R[2,2]) == 1): # R[0,2] == 0 and R[1,2] == 0
        if(R[2,2]==1): #theta=0
            theta = 0
            phi = 0
            psi = np.arctan2(-R[0,1],R[1,1])
            #psi = np.arctan2(R[1,0],R[0,0]) # otra forma
        elif(R[2,2]==1):
            theta = pi
            phi = np.arctan2(-R[0,1],R[1,1])
            #phi = np.arctan2(-R[1,0],-R[0,0]) # otra forma
            psi = 0
    else:
        print("ERROR")

    return phi, theta, psi

def matriz_R_euler(phi, theta, psi):
    
    s_ph = sin(phi)
    c_ph = cos(phi)
    s_th = sin(theta)
    c_th = cos(theta)
    s_ps = sin(psi)
    c_ps = cos(psi)
    
    # VERIFICAMOS EL ANGULO "phi"
    if(abs(phi)==pi or abs(phi) == 2*pi):
        s_ph = 0
    elif(abs(phi) == pi/2):
        c_ph = 0

    # VERIFICAMOS EL ANGULO "theta"
    if(abs(theta)==pi or abs(theta) == 2*pi):
        s_th = 0
    elif(abs(theta) == pi/2):
        c_th = 0
        
    # VERIFICAMOS EL ANGULO "psi"
    if(abs(psi)==pi or abs(psi) == 2*pi):
        s_ps = 0
    elif(abs(psi) == pi/2):
        c_ps = 0
    
    r = [[c_ph*c_th*c_ps - s_ph*s_ps, -c_ph*c_th*s_ps - s_ph*c_ps, c_ph*s_th],
         [s_ph*c_th*c_ps + c_ph*s_ps, -s_ph*c_th*s_ps + c_ph*c_ps, s_ph*s_th],
         [-s_th*c_ps, s_th*s_ps, c_th]]
    R = np.array(r)
    
    return R

def matriz_R_DH(theta, alpha):
    
    s_th = sin(theta)
    c_th = cos(theta)
    s_al = sin(alpha)
    c_al = cos(alpha)
    
    # VERIFICAMOS EL ANGULO "theta"
    if(abs(theta)==pi or abs(theta) == 2*pi):
        s_th = 0
    elif(abs(theta) == pi/2):
        c_th = 0
        
    # VERIFICAMOS EL ANGULO "alpha"
    if(abs(alpha)==pi or abs(alpha) == 2*pi):
        s_al = 0
    elif(abs(alpha) == pi/2):
        c_al = 0
    
    r = [[c_th, -s_th*c_al, s_th*s_al],
         [s_th, c_th*c_al, -c_th*s_al],
         [0, s_al, c_al]]
    R = np.array(r)
    
    return R

def cinematica_directa(q):  
    ## Parametros D-H
    theta = [q[0], q[1], q[2], q[3], q[4], q[5], q[6]]
    d = [0.36, 0, 0.42, 0, 0.4, 0, 0.126]
    a = [0, 0, 0, 0, 0, 0, 0]
    alpha = [-pi/2, pi/2, -pi/2, pi/2, -pi/2, pi/2, -pi/2]

    A1 = matriz_homogenea_DH(q[0], d[0], a[0], alpha[0])
    A2 = matriz_homogenea_DH(q[1], d[1], a[1], alpha[1])
    A3 = matriz_homogenea_DH(q[2], d[2], a[2], alpha[2])
    A4 = matriz_homogenea_DH(q[3], d[3], a[3], alpha[3])
    A5 = matriz_homogenea_DH(q[4], d[4], a[4], alpha[4])
    A6 = matriz_homogenea_DH(q[5], d[5], a[5], alpha[5])
    A7 = matriz_homogenea_DH(q[6], d[6], a[6], alpha[6])

    T1 = A1
    T2 = np.dot(T1,A2)
    T3 = np.dot(T2,A3)
    T4 = np.dot(T3,A4)
    T5 = np.dot(T4,A5)
    T6 = np.dot(T5,A6)
    T7 = np.dot(T6,A7)
    H = T7
    
    pos_inicial = np.reshape(np.array([0, 0, 0, 1]), (4,1))
    pos_final = np.dot(H,pos_inicial)
    x = float(pos_final[0])
    y = float(pos_final[1])
    z = float(pos_final[2])
    phi, theta, psi = angulos_euler(H)
    
    return H, x, y, z, phi, theta, psi

def error_rms(pos_deseada, orientacion, q):
    
    H, x, y, z, phi, theta, psi = cinematica_directa(q)
    error_x = pos_deseada[0]-x
    error_y = pos_deseada[1]-y
    error_z = pos_deseada[2]-z
    error_rms = float(np.sqrt((error_x**2 + error_y**2 + error_z**2)/3))

    return error_rms

def flag_cinematica_inversa(pos_deseada, orientacion, flag1, flag2, flag4, flag6):
    
    # Parametros D-H
    #theta = [q1, q2, q3, q4, q5, q6, q7]
    q = [0, 0, 0, 0, 0, 0, 0]
    d = [0.36, 0, 0.42, 0, 0.4, 0, 0.126]
    a = [0, 0, 0, 0, 0, 0, 0]
    alpha = [-pi/2, pi/2, -pi/2, pi/2, -pi/2, pi/2, -pi/2]

    R = matriz_R_euler(orientacion[0], orientacion[1], orientacion[2])

    ## Se define un angulo para q7 = q[6] inicial
    q[6] = 0 # q7=0

    ## Se halla la junta q4 = q[3]
    R67 = matriz_R_DH(q[6], alpha[6])
    R06 = np.dot(R, np.linalg.inv(R67))
    P07 = np.reshape(np.array(pos_deseada), (3,1))
    P67 = np.reshape(np.array([0, 0, d[6]]), (3,1))
    d06 = P07 - np.dot(R06, P67)
    L = np.sqrt((d06[0])**2 + (d06[1])**2 + (d06[2] - d[0])**2)
    cq4 = (L**2 - d[2]**2 - d[4]**2)/(2*d[2]*d[4])

    q4 = np.arctan2(flag4*np.sqrt(max(1-cq4**2,0)),cq4)   # q4 -> [0,pi]
    #q4 = np.arctan2(-np.sqrt(max(1-cq4**2,0)),cq4)  # q4 -> [-pi,0]
    q[3] = float(q4) # q4 -> [0,pi]

    ## Se define un angulo para q3 = q[2]
    q[2] = 0

    ## Se halla la junta q1 = q[0]
    sk1_N = d[4]*sin(q[2])*sin(q[3])
    sk1_D = np.sqrt((d06[0])**2+(d06[1])**2)
    q1 = np.arctan2(d06[1],d06[0])-np.arctan2(sk1_N,flag1*np.sqrt(max(sk1_D**2-sk1_N**2,0))) # q1 -> [-pi/2,pi/2]
    #q1 = np.arctan2(d06[1],d06[0])-np.arctan2(sk1_N,-np.sqrt(max(sk1_D**2-sk1_N**2,0))) # q1 -> [-pi, -pi/2]U[pi/2, pi]
    q[0] = float(q1)

    ## Se halla la junta q2 = q[1]
    k1 = d[4]*cos(q[2])*sin(q[3])
    k2 = d[4]*cos(q[3]) + d[2]
    k3 = d06[2] - d[0]
    ck2 = k3/np.sqrt(k1**2 + k2**2)

    q2 = np.arctan2(flag2*np.sqrt(max(1-ck2**2,0)), ck2) - np.arctan2(k1,k2)  # q2 -> [0,pi]
    #q2 = np.arctan2(-np.sqrt(max(1-ck2**2)), ck2) - np.arctan2(k1,k2)  # q2 -> [-pi,0]
    q[1] = float(q2) 

    ## Se halla la junta q6 = q[5], q5 = q[4] y q7 = q[6]
    A1 = matriz_homogenea_DH(q[0], d[0], a[0], alpha[0])
    A2 = matriz_homogenea_DH(q[1], d[1], a[1], alpha[1])
    A3 = matriz_homogenea_DH(q[2], d[2], a[2], alpha[2])
    A4 = matriz_homogenea_DH(q[3], d[3], a[3], alpha[3])
    T1 = A1
    T2 = np.dot(T1,A2)
    T3 = np.dot(T2,A3)
    T4 = np.dot(T3,A4)

    R4 = T4[0:3,0:3]
    R47 = np.dot(np.linalg.inv(R4), R)
    M = R47

    ck6 = -M[2,1]
    q[5] = np.arctan2(flag6*np.sqrt(max(1-ck6**2,0)), ck6) # q6 -> [0, pi]
    #q[5] = np.arctan2(-np.sqrt(max(1-ck6**2,0)), ck6) # q6 -> [-pi, 0]
    
    if sin(q[6])>0:
        q[4] = np.arctan2(-M[1,1], -M[0,1]) # q5
        q[6] = np.arctan2(M[2,2], -M[2,0]) # q7
    else:
        q[4] = np.arctan2(M[1,1], M[0,1]) # q5
        q[6] = np.arctan2(-M[2,2], M[2,0]) # q7

    return q

def cinematica_inversa(pos_deseada, orientacion, error=0.01):
    
    for a in [1,-1]:
        for b in [1,-1]:
            for c in [1,-1]:
                for d in [1,-1]:
                    q1 = flag_cinematica_inversa(pos_deseada, orientacion, a, b, c, d)
                    if (error_rms(pos_deseada, orientacion, q1) < error): break
                if (error_rms(pos_deseada, orientacion, q1) < error): break
            if (error_rms(pos_deseada, orientacion, q1) < error): break
        if (error_rms(pos_deseada, orientacion, q1) < error): break

    #print("{} {} {} {}".format(a,b,c,d))  
    q = q1
    
    return q

