import os
import numpy as np

wave=0
flux=0
sigup=0
normwave=0
normflux=0
normsig=0
fitidx=0
fitcoeff=0
fitcovmtx=0
fiterrors=0
origpars=0
origparinfo=0
filename=''
field=''
homedir=os.path.expanduser('~')
spectrum = []

Todd = True
if Todd:
    lsf='COS_LP1'
    instr=['COS','COS','COS','COS','STIS']
    gratings=['G130M','G160M','G185M','G225M','E230M']
    slits=['NA','NA','NA','NA','0.2x0.2']
    lsfranges=np.array([[1100,1400],[1400,1800],[1800,2100],[2100,2278],[2278,4100]])
    lps=['1','1','1','1','1']
    cen_wave=['1327','1600','1953','2250','2707']

else:  # this is for the casual user
    lsf='COS_LP1'
    instr=['COS','COS']
    gratings=['G130M','G160M']
    slits=['NA','NA']
    lsfranges=np.array([[1100,1450],[1450,1800]])
    lps=['1','1']
    cen_wave=['1327','1600']


echarge=4.803204505713468e-10
m_e=9.10938291e-28
c=29979245800.0

lyseries=np.array([913.8260,914.0390,914.2860,914.5760,914.919,915.329,915.824,916.429,917.1806,918.1294,919.3514,920.9631,923.1504,926.2257,930.7483,937.8035,949.7431,972.5368,1025.7223,1215.6701])
NIseries=[1199.5496, 1200.2233, 1200.7098]
NIIseries=[533.5099, 644.6337, 671.4102, 915.6131, 1083.9937]
NIIIseries=[374.204,684.996,685.513,763.34,989.799]
NIVseries=[247.205, 765.148, 1486.496]
NVseries=[1238.821, 1242.804]
OIseries=[791.5136, 791.9732,877.7983, 877.8787,971.7371, 971.7376, 971.7382,988.5778,988.6549, 988.7734,1302.1685]
OIIseries=[376.693, 376.745, 387.651, 387.901, 388.058, 391.912, 391.943, 392.002, 418.598, 418.881, 419.065, 429.918, 430.041, 430.177, 539.0855, 539.5489, 539.8544, 832.7572, 833.3294, 834.4655]
OIIIseries=[228.825, 248.538, 255.155, 262.662, 263.692, 264.257, 266.967, 275.281, 280.234, 303.411, 305.596, 374.005, 475.153, 507.391, 702.332, 832.927]
OIVseries=[238.36, 279.631, 553.33, 554.075, 608.398, 787.711, 1397.232, 1399.78]
OVseries=[629.73, 1218.344]
OVIseries=[1031.9261, 1037.6167]
SiIIseries=[889.7228, 899.4063, 989.8731, 1020.6989, 1190.4158, 1193.2897, 1260.4221, 1304.3702, 1526.7066, 1808.0126]
SiIIIseries=[408.378, 423.817, 426.644, 437.255, 466.129, 566.613, 1206.5, 1892.03]
SiIVseries=[300.517, 300.531, 310.234, 310.257, 327.137, 327.181, 361.56, 361.659, 457.818, 458.155, 1393.755, 1402.77]
SIIseries=[558.755, 558.924, 559.131, 571.156, 571.364, 571.779, 574.397, 576.057, 576.978, 593.507, 593.835, 594.475, 600.661, 602.446, 603.43, 640.412, 640.902, 641.767, 662.267, 664.315, 665.519, 763.657, 764.42, 765.693, 906.885, 910.484, 912.735, 1250.578, 1253.805, 1259.518]
SIIIseries=[476.351, 480.533, 484.172, 677.746, 681.47, 698.731, 724.289, 1012.501, 1190.191]
SIVseries=[368.99, 391.558, 551.17, 657.34, 744.907, 748.4, 809.668, 1062.662, 1398.05, 1404.8]
SVseries=[286.094, 786.48, 1199.134]
SVIseries=[248.99, 249.27, 933.378, 944.523]
CIseries=[945.191, 1277.2454, 1280.1353, 1328.8333, 1560.3092, 1613.3763, 1656.9283]
CIIseries=[438.773, 438.824, 466.352, 466.407, 530.274, 538.406, 543.257, 549.3195, 549.3785, 551.68, 560.2394, 576.8748, 594.8, 635.9945, 687.0526, 858.0918, 903.6235, 903.9616, 1036.3367, 1334.5323, 2324.214, 2325.403]
CIIIseries=[270.324, 274.051, 280.043, 288.423, 291.3261, 310.1697, 322.5741, 386.2028, 977.02, 1908.734]
CIVseries=[244.907, 244.911, 312.422, 312.453, 1548.195, 1550.77]
FeIIseries=[2260.7805, 2344.214, 2367.5905, 2374.4612, 2382.765, 2586.65, 2600.1729]
FeIIIseries=[734.296, 735.247, 735.344, 808.078, 813.379, 814.137, 823.256, 824.799, 826.387, 844.288, 850.907, 854.2, 858.551, 858.609, 859.723, 1122.524, 1207.05, 1214.562]
FeIVseries=[525.689, 526.293, 526.634]
FeVseries=[363.444, 365.44, 386.261, 387.984]
FeXIVseries=[252.19, 257.385, 274.203, 334.178]
FeXVIseries=[335.407, 360.798]
NeIIseries=[445.0397, 446.2556, 460.7284]
NeIIIseries=[313.059, 345.448, 488.108, 489.505]
NeIVseries=[334.056, 541.127, 542.073, 543.891]
NeVseries=[329.151, 357.95, 480.41, 568.42]
NeVIseries=[399.82, 401.14, 433.176, 558.59, 993.0, 997.4]
NeVIIseries=[465.221, 895.18]
NeVIIIseries=[770.409, 780.324]
MgIXseries=[368.071, 705.72]
MgXseries=[609.79, 624.95]
ArIVseries=[451.2, 451.87, 452.91, 840.03, 843.77, 850.6]
HeIseries=[506.2, 506.5702, 507.0576, 507.7178, 508.6431, 509.9979, 512.0982, 515.6165, 522.2128, 537.0296, 584.334]
HeIIseries=[229.431, 229.736, 230.139, 230.686, 231.454, 232.584, 234.347, 237.331, 243.027, 256.317, 303.7822]

multiplets=[lyseries,NIseries,NIIseries,NIIIseries,NIVseries,NVseries,OIseries,OIIseries,OIIIseries,OIVseries,OVseries,OVIseries,FeIIseries,FeIIIseries,FeIVseries,FeVseries,FeXIVseries,FeXVIseries,NeIIseries,NeIIIseries,NeIVseries,NeVseries,NeVIseries,NeVIIseries,NeVIIIseries,MgIXseries,MgXseries,ArIVseries,HeIseries,HeIIseries,CIseries,CIIseries,CIIseries,CIVseries,]

wavegroups=[]
wgidxs=[]
uqwgidxs=[]
lsfs=[]
fgs=[]

#Atomic data
lams=np.array([])
fosc=np.array([])
gam=np.array([])

outputdir='./'

VPparoutfile=outputdir+field+'_VP.dat'
VPmodeloutfile=outputdir+field+'VPmodel.fits'
contoutfile=outputdir+'continua.dat'
largeVPparfile=outputdir+'_VP_log.dat'
defaultcol=13.1
defaultb=20.0
defaultvlim=1000.
lowblim=4.
upperblim=85.
upperblim_HI=210.

# Plotting
ylim = (-0.2, 1.4)
xtick_fontsize = 'small'
ytick_fontsize = 'small'
xy_fontsize = 'small'
x_labelpad = 0
y_labelpad = -3
label_ypos = 0.2 * (ylim[0]+ylim[1])
label_fontsize = 8
