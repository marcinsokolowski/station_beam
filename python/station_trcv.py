#!/usr/bin/python
from __future__ import print_function
import sys
import math
import numpy
from scipy.interpolate import interp1d

from optparse import OptionParser

def db2num( val_db ) :
   val_num = 10.00**( val_db / 10.0 )
   
   return val_num

# Based on values from Budi : /home/msok/Desktop/EDA/data/Budi/EDA_LNA_comparison/EDA_LNA_comparison_Marcin_2017_02_22.csv
#      awk '{printf("%s,",$1);}' EDA_LNA_comparison_Marcin_2017_02_22.txt
#      
def eda2_lna_gain_db( freq_mhz ) :
   # 49 was added so that results for 50 MHz make sense :
   freq_arr = [49, 50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650]
   lna_gain_db = [9.605079504, 9.605079504,11.03451552,12.00310915,12.87530424,13.66458303,14.59206758,15.20004633,15.81631652,16.57700878,17.04504641,17.4530421,17.94119443,18.29400522,18.64742568,18.99786667,19.1927476,19.33077469,19.55665415,19.82231412,20.05328333,20.01965598,20.25221114,20.40278841,20.43272563,20.48501645,20.57134669,20.6484693,20.72924232,20.69060146,20.79451705,20.72767332,20.78805191,20.88447966,20.81845264,20.84605345,20.96364294,21.02841924,21.03875741,20.86557661,21.01645346,21.00397445,21.00794548,20.96025789,21.06591976,21.06959974,20.97894178,21.03926669,21.18819212,21.15661213,21.10576505,21.04062147,21.02049276,21.13065673,21.18310726,21.1494003,21.23756747,21.15425812,21.27121927,21.23367783,21.18590315,21.25894461,21.35938183,21.24941471,21.2457331,21.13270465,21.05061006,21.2878487,21.26380057,21.30918524,21.22981428,21.27843001,21.25678546,21.13101944,21.34725711,21.26598385,21.24300716,21.25646715,21.19995247,21.31016549,21.36595784,21.31436298,21.23127767,21.24182649,21.22408933,21.18858055,21.18252406,21.31899335,21.29157879,21.26990699,21.23491372,21.17830177,21.2312951,21.15408903,21.21865198,21.18213178,21.11416074,21.10051057,21.13797436,20.99505689,20.99465571,21.01001425,21.07220168,20.91550143,20.84248206,20.90911656,20.79285463,20.82641257,20.77193833,20.60304284,20.70168979,20.7923329,20.74332209,20.71229048,20.71419611,20.57517084,20.53713666,20.5910783,20.60171482,20.38349074,20.42532898,20.35405356,20.38371749,20.36302448,20.31521109,20.21948596,20.13346744,20.31200363,20.12268727,20.07306271,20.09158655,19.9740745,19.97497133,19.8918078,19.91257659,19.94306702,19.78789876,19.83310032,19.78268469,19.64582653,19.74286006,19.64862211,19.50780087,19.58466507,19.53413535,19.51661912,19.36417446,19.35174632,19.37785439,19.33774608,19.22187731,19.24780933,19.30682268,19.1371701,19.18335652,19.16838572,19.02875004,19.04273242,18.95725781,18.98774777,18.8674208,19.00153798,18.84735881,18.86215619,18.80189323,18.74442182,18.74756168,18.74632534,18.63803679,18.70541864,18.70671949,18.63756359,18.59120836,18.55837469,18.49926607,18.50797131,18.52480355,18.43462456,18.4294043,18.37367211,18.28010532,18.251653,18.24348731,18.1778639,18.26330685,18.16781137,18.14388727,18.25384816,18.05312782,18.12694256,18.04981574,18.02437818,18.07282708,18.1444027,18.05991263,18.00279867,18.01583099,17.98916911,17.92513427,17.84324583,17.92516187,17.91812968,17.8920843,17.81488047,17.74206942,17.79835054,17.86894264,17.74069144,17.70476884,17.76796544,17.72365227,17.71328266,17.77821747,17.81615389,17.62618148,17.6468922,17.58101943,17.71626977,17.66875583,17.71639647,17.65727591,17.64349273,17.63530495,17.70615123,17.60013171,17.60024122,17.68090062,17.71512435,17.63662816,17.63690703,17.68543104,17.64885749,17.5502575,17.65905747,17.65262169,17.612775,17.58670317,17.68127316,17.68947478,17.73594805,17.72298306,17.82434602,17.6183226,17.79222567,17.75450593,17.72091597,17.75896537,17.80798817,17.73048324,17.89163527,17.73154054,17.73695321,17.87860697,17.7899376,17.89676563,17.94776504,17.84097877,17.87279891,17.97810731,17.93045113,17.93477882,17.97810386,17.92704977,17.97460896,18.04433376,18.09445622,17.95757422,18.11389069,18.06881443,18.08145777,18.23010367,18.09855644,18.08460046,18.24036162,18.23494619,18.26321571,18.18595116,18.34624281,18.36132165,18.38751,18.287058,18.35693174,18.32024521,18.3372692,18.54168736,18.50734984,18.47346192,18.49828763,18.54180785,18.6248263,18.65329682,18.66730654,18.70704846,18.7444835,18.70099517,18.81259454,18.92669263,18.89860224,19.09362459,18.77869359,18.7777383,18.82928545,18.90425779,19.04544793,19.0473181,19.02879157,19.16788356,19.21911325,19.07745657,19.25826503,19.23550131,19.29683114,19.27589453,19.43331668,19.29782679,19.32901814,19.41098107,19.53416138,19.60735117,19.55424815,19.49281101,19.63977841,19.67078198,19.69516313,19.7461468,19.79074561,19.78972399,19.81228526,19.81955658,19.96404773,19.87885777,19.95246346,19.93666889,19.94161093,20.05872015,20.12710505,20.10996774,20.09944281,20.20504613,20.17460593,20.21880734,20.3510194,20.20599369,20.31260715,20.15107945,20.34090427,20.38659683,20.39884135,20.32102356,20.3811628,20.39234166,20.44055156,20.38447773,20.43609748,20.54390459,20.44029915,20.45744847,20.48357014,20.42490441,20.50300829,20.41029967,20.54323908,20.58954463,20.42616971,20.45698106,20.47686166,20.49328035,20.56950041,20.45248547,20.51620407,20.36521404,20.46763153,20.44264701,20.40747878,20.38795384,20.40455263,20.41822692,20.34185965,20.33505076,20.32517606,20.17405417,20.26908089,20.2066725,20.20383006,20.193328,20.15711882,20.10126536,20.06962938,20.02342236,20.03130727,20.04339041,19.83004028,19.84734956,19.8230425,19.82363804,19.80213875,19.7660997,19.62819711,19.50156698,19.54741814,19.40418184,19.5315688,19.4558576,19.38311988,19.22626574,19.1757151,19.21435866,19.19905422,19.06308453,19.01257725,19.03100578,18.89384554,18.97641616,18.80424198,18.87821146,18.74707601,18.59272058,18.6289792,18.56796329,18.56745448,18.50577073,18.52101985,18.39590434,18.35637683,18.21644862,18.26788588,18.23197965,18.051489,18.0626333,18.00907116,18.12388702,17.85101587,17.97224087,17.84346985,17.75318774,17.70741646,17.73465671,17.67098987,17.63084555,17.53200395,17.44242259,17.41544296,17.35771924,17.41270032,17.27093035,17.25743394,17.13626385,17.22464822,17.09216549,17.07660044,16.8521157,16.9967049,16.90783091,16.96213896,16.83775207,16.88138649,16.67559739,16.62178529,16.63296946,16.53175421,16.48462581,16.53973612,16.4297137,16.3482631,16.33958727,16.30497681,16.24823316,16.21489363,16.17405492,16.09015195,16.11803083,15.95864686,16.03753419,15.86460183,15.82943409,15.75316556,15.80816924,15.7409121,15.68171871,15.70085424,15.69307462,15.58151465,15.68344449,15.52179622,15.6558834,15.43457798,15.43684035,15.50839516,15.56005626,15.38896236,15.35868547,15.49164344,15.24407692,15.33134968,15.27798492,15.25716566,15.19828627,15.08024949,14.95174184,14.85883818,14.89582439,14.81194749,14.7195624,14.69437054,14.69503321,14.7535969,14.65146094,14.57623544,14.72752293,14.60542195,14.57806476,14.65132141,14.60387094,14.56268799,14.50276442,14.54098724,14.42762028,14.54815744,14.46416894,14.42676494,14.41052527,14.43693176,14.32499456,14.27346364,14.37957709,14.22951782,14.29193628,14.27875712,14.28279964,14.20901835,14.21540508,14.1571718,14.20454074,14.09116181,14.18230469,14.0407206,14.0165534,14.05426597,13.95409985,14.04579286,13.87890277,14.03013758,13.92112221,13.8416521,13.9704673,13.89797864,13.77953873,13.83398131,13.85850691,13.76421088,13.83676199,13.87036358,13.75878235,13.9000158,13.82128952,13.73015675,13.8688501,13.72128559,13.85416925,13.70103055,13.7431367,13.72323175,13.78536671,13.56380074,13.73344337,13.69967724,13.64117484,13.63096654,13.67406216,13.79450907,13.61927781,13.72881562,13.69614091,13.61236544,13.57471721,13.63158441,13.57887348,13.7367574,13.69813475,13.54577761,13.59747822,13.61661789,13.59816776,13.6448558,13.41620817,13.55390495,13.63219641,13.5873037,13.64188166,13.53219296,13.60532578,13.61689383,13.55374843,13.53528873,13.70905433,13.56987374,13.6115553,13.53619308,13.67149637,13.56016505,13.51127688,13.4892625]
   
   lna_gain_db_cubic = interp1d( freq_arr, lna_gain_db, kind='cubic')

   ret_val = float( lna_gain_db_cubic( freq_mhz ) )
   print("DEBUG : EDA2 LNA GAIN at %.2f MHz = %.4f [dB]" % (freq_mhz,ret_val))
   return ret_val

# see ~/Desktop/EDA/loogbook/Cath/Budi_trcv and ~/Desktop/EDA/loogbook/Cath/eda_aeff_cath.odt
# Budi's e-mail on 2016-12-20 
def trcv_budi( freq_mhz ):
   t_lna =  +724.02940227341059653554111719131469726562500000000000 -18.91813674541278089691331842914223670959472656250000*freq_mhz \
             +0.23407541076266874524591798945039045065641403198242*freq_mhz*freq_mhz -0.00233671195682318014533174199698351003462448716164*freq_mhz*freq_mhz*freq_mhz \
             +0.00002254381274811321963493990216953477556671714410*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000016409668824170788921994429236661527937712890*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000076016725276938624039072744188781996510507*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz - 0.00000000000210443395476647521171394324012829569895*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000317128431524229797923526252249329818*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz - 0.00000000000000000179318498778901811500861463716146*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000034683893039093007945200578004*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             - 0.00000000000000000000000036757731954965862822568812*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000262855743668854032145143 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000001701169748990267910934 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000001832975188922028666 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000010411258968305159 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000036864325213161 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000000061278189979 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000000000000000000000000000000147452601 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000000000000000000000000000000000001736900 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   print("INFO : using Budi's data for T_rcv ( %.2f MHz) = %.2f" % (freq_mhz,t_lna))

   return t_lna   

# fitted from POWER vs. MODEL data in LST range 18-20 hours only see Ill. Illustration 6 in eda_drift_scan_data_201612.odt for details:
def trcv_fit_data_vs_mode( freq_mhz ):
   t_lna = +40641.44874588610400678589940071105957031250000000000000  -1151.44119960312536932178772985935211181640625000000000*freq_mhz \
           +10.28736969420266333941071934532374143600463867187500 * freq_mhz*freq_mhz -0.00920101512017338675486488597243805998004972934723 * freq_mhz*freq_mhz*freq_mhz -0.00024102523369100168010650819816476086998591199517 * freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           -0.00000019651325627217189202790260330738769667391352 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000820101610522132312859120128458942367011275 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           +0.00000000001047686030411538916735193183367932131915 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000000024494646525009589774756356651373232348 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
           +0.00000000000000052438247226051691783257993415371231 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   return t_lna

# Fitted - see :
# cd /home/msok/Desktop/EDA/data/2016-07/20160703
# root [4] .x plot_tlna_MHz.C("trcv_vs_freq.selected",NULL,0,5000)
# t_lna =  +36853.38919790466752601787447929382324218750000000000000 * 1 -1377.64458069924739902489818632602691650390625000000000 * 1*freq_mhz +20.37196908756931890138730523176491260528564453125000 * 1*freq_mhz*freq_mhz -0.13515251806813355361924777753301896154880523681641 * 1*freq_mhz*freq_mhz*freq_mhz +0.00022215700342832035230414278181854115246096625924 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000169842384409759492976880849990362065682347747 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000464297999395697167928294570157440879221156 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000003336389317067374617426693716906293402102 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000000017598078793649799250229898970650520328 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000000000023477806740472191603243644100388291 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz;
# Cubic interpolation of data points see  T_rcv from lightcurve fits see Ill.1 and Ill.2 in trcv_from_20160703_drift_scan.odt for details
def trcv_fit_lightcurve_20160703_cubic( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' trcv_vs_freq.selected
   # awk '{printf("%.2f,",$2);}' trcv_vs_freq.selected
   x=[55.00,65.00,75.00,85.00,95.00,105.00,115.00,145.00,155.00,165.00,175.00,185.00,195.00,205.00,215.00,225.00,235.00]
   y=[3831.09,1744.61,1022.83,695.13,489.77,340.59,232.78,83.88,59.21,45.78,35.69,30.60,27.26,30.52,27.69,28.39,93.70]
      
   tlna_cubic = interp1d(x, y, kind='cubic')
   
   return float( tlna_cubic(freq_mhz) )

# EDA/loogbook/haslam_vs_angelica.odt Ill. 24 red curve
# /home/msok/Desktop/EDA/data/2016-07/20160703/ANGELICA/trcv_vs_freq.selected
def trcv_angelica_data_vs_time( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' ANGELICA/trcv_vs_freq.selected   
   # awk '{printf("%.2f,",$2);}' ANGELICA/trcv_vs_freq.selected
   x=[55.00,65.00,75.00,85.00,95.00,105.00,125.00,155.00,175.00,205.00,225.00]
   y=[3322.52,1608.02,1007.08,824.80,638.77,232.16,245.80,68.14,47.27,33.88,27.79]
   
   tlna_cubic = interp1d(x, y, kind='cubic')
   
   return float( tlna_cubic(freq_mhz) )
   
def trcv_angelica_data_vs_time_powerlawfit( freq_mhz ) :
   t_rcv = 384.3 * math.pow((100.00/freq_mhz),3.571)
   
   if freq_mhz >= 160 and False :      
      t_rcv = trcv_budi(freq_mhz)/corr_factor_adrian(freq_mhz)
      
   return t_rcv   
      

# EDA/loogbook/haslam_vs_angelica.odt Ill. 24 green curve
# /home/msok/Desktop/EDA/data/2016-07/20160703/HASLAM/trcv_vs_freq.data_vs_model_all
def trcv_haslam_data_vs_model( freq_mhz ):
   # cd ~/Desktop/EDA/data/2016-07/20160703
   # awk '{printf("%.2f,",$1);}' HASLAM/trcv_vs_freq.data_vs_model_all
   # awk '{printf("%.2f,",$2);}' HASLAM/trcv_vs_freq.data_vs_model_all
   x=[55.00,65.00,75.00,85.00,95.00,105.00,115.00,145.00,155.00,165.00,175.00,185.00,195.00,205.00,215.00,225.00,235.00]
   y=[3849.12,1744.54,1013.68,682.12,473.69,327.36,222.33,81.70,57.44,45.02,35.46,30.50,27.35,31.30,28.75,29.14,45.67]
   
   tlna_cubic = interp1d(x, y, kind='cubic')
   return float( tlna_cubic(freq_mhz) )

# as in loogbook/paper/eda_lightcurve_and_trcv.odt Ill. 13:
# /home/msok/Desktop/EDA/data/2016-07/20160703/HASLAM/WithSun_and_BeamformingErrors/BIGHORNS/images/final/paper/trcv_vs_freq_lst_13-17_hours.png
# root [1] .x plotspec_freq_errxy.C("trcv_vs_freq_lst_range_13_17_hours.txt")
# Fitted : [0]*(150/x)^[1]
def trcv_from_skymodel_with_err( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too 
   # not to go below 50K (as fitting might be tricky)
#   if trcv < 50 or freq_mhz>160 :
#   if freq_mhz > 160 :
#      return 50

   print("INFO : using trcv_from_skymodel_with_err")

   if use_cubic :
      x=[          55.0,          65.0,        75.0,         85.0,         95.0,       105.0,        115.0,        125.0,        135.0,        145.0,       155.0,       165.0,       175.0,       185.0,       195.0,       205.0]
      y=[ 3722.55055596, 1421.83480972, 655.8033454, 409.69506365, 356.59552554, 293.4326318, 200.26593869, 163.91940732, 157.46121101, 112.19530177, 51.14849704, 68.56693587, 28.33318976, 43.72596282, 34.41633474, 62.76422256]   
      tlna_cubic = interp1d(x, y, kind='cubic')
      trcv = float( tlna_cubic(freq_mhz) )
   else :
      index=3.46114
      trcv = 82.4708*math.pow( (150.00/freq_mhz), index )
   
   if trcv < 50 or freq_mhz>160 :
      trcv = 50
      
   # see /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity.odt or /home/msok/Desktop/EDA/doc/T_rcv   
#   if freq_mhz >= 240 :
#      trcv = 70 

   # see /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity.odt and /home/msok/Desktop/EDA2/logbook/20200108_expected_sensitivity_FINAL.odt
   if freq_mhz > 160 :
      trcv = 30.2000 + freq_mhz*0.100125
   
   return trcv


# AAVS2 : 
# see : /home/msok/Desktop/AAVS2/logbook/20200128_receiver_temperature.odt
def trcv_aavs2( freq_mhz , use_cubic=True, za_deg=None ): # default was power law fit but I wanted to be able to use cubic fit too 
   # not to go below 50K (as fitting might be tricky)
#   if trcv < 50 or freq_mhz>160 :
#   if freq_mhz > 160 :
#      return 50


   # at the moment return the same as for EDA-1 (or EDA-2)
#   return trcv_from_skymodel_with_err( freq_mhz , use_cubic )

   trcv = 33.85 + 0.0148*freq_mhz 
   
   using_cubic = False
   if freq_mhz < 100.00 :
      x=[10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,42.0,44.0,46.0,48.0,50.0,52.0,54.0,56.0,58.0,60.0,62.0,64.0,66.0,68.0,70.0,72.0,74.0,76.0,78.0,80.0,82.0,84.0,86.0,88.0,90.0,92.0,94.0,96.0,98.0,100.0,102.0,104.0,106.0,108.0,110.0,112.0,114.0,116.0,118.0,120.0,122.0,124.0,126.0,128.0,130.0,132.0,134.0,136.0,138.0,140.0,142.0,144.0,146.0,148.0,150.0,152.0,154.0,156.0,158.0,160.0,162.0,164.0,166.0,168.0,170.0,172.0,174.0,176.0,178.0,180.0,182.0,184.0,186.0,188.0,190.0,192.0,194.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,292.0,294.0,296.0,298.0,300.0,302.0,304.0,306.0,308.0,310.0,312.0,314.0,316.0,318.0,320.0,322.0,324.0,326.0,328.0,330.0,332.0,334.0,336.0,338.0,340.0,342.0,344.0,346.0,348.0,350.0,352.0,354.0,356.0,358.0,360.0,362.0,364.0,366.0,368.0,370.0,372.0,374.0,376.0,378.0,380.0,382.0,384.0,386.0,388.0,390.0,392.0,394.0,396.0,398.0,400.0,402.0,404.0,406.0,408.0,410.0]
      y=[7107100.0000,654850.0000,206830.0000,376900.0000,361870.0000,361320.0000,384580.0000,408710.0000,370030.0000,257710.0000,194070.0000,108050.0000,55808.0000,29746.0000,17932.0000,11791.0000,7751.7000,5199.5000,3493.5000,2188.8000,1356.3000,747.7300,401.0200,210.2100,118.1100,82.2620,72.2770,68.3100,64.0260,57.4240,49.0060,43.9750,40.7470,41.5500,44.0500,46.5410,48.0550,48.9760,46.7510,42.1360,39.1750,33.6010,32.7680,32.4490,32.4100,35.3490,33.8180,35.4250,36.2200,37.5340,36.2480,36.4220,35.7270,35.0100,35.8200,35.4200,37.2410,37.0950,36.3850,35.0540,33.2270,33.0210,35.3170,36.1650,38.3440,38.6990,36.5530,35.4320,34.9590,34.6170,36.0870,37.4090,38.1410,39.4630,37.1560,36.1880,36.3860,35.1380,35.4160,35.5970,37.9060,38.6590,38.3350,36.6420,37.6530,36.6690,34.2690,34.3780,33.7400,35.3170,36.9150,37.8910,37.6220,36.9020,37.7700,40.9420,35.1330,35.2780,33.5290,32.9030,34.3800,34.4300,37.6530,39.4540,39.4440,41.0410,40.7670,40.4160,40.1550,39.6900,38.7630,36.3620,35.3780,34.3560,33.9780,33.9440,33.6700,32.6710,34.1480,34.1550,34.2700,35.1280,34.9730,34.7650,34.0960,35.0130,35.2690,36.0910,36.1660,38.7710,38.3090,39.7370,39.9530,41.9870,42.0490,37.8430,40.1820,42.6590,40.5190,39.8170,40.5610,38.7800,37.9280,37.0910,37.1220,36.7460,38.1390,37.4310,37.3780,39.4030,39.4210,39.3430,40.6120,41.0690,40.2540,36.9340,40.3270,39.4150,39.3720,39.0670,38.4750,36.8540,36.4600,37.2850,37.3790,37.4560,38.3510,36.9400,38.7410,39.1510,38.8910,42.3500,41.2200,41.3130,40.9190,37.1090,39.4140,40.1410,38.1900,37.2650,37.8390,37.3610,36.5130,36.7570,37.5530,37.1360,39.0620,39.9190,42.8730,44.6690,45.2820,45.9850,47.6020,47.5390,49.5910,47.3790,49.9540,49.2800,49.2160,49.2670,47.3600]
      
      trcv_interpol = interp1d(x, y, kind='cubic')
      trcv = float( trcv_interpol( freq_mhz ) )
      using_cubic = True

   if za_deg is not None :
      print("za_deg is not None:")
      trcv_vs_za_fit_params={}
      trcv_vs_za_fit_params[80]  = numpy.array( [ 44.37610000 , 0.01785778 , 0.00073852 ] )
      trcv_vs_za_fit_params[110] = numpy.array( [ 39.43860000 , 0.02064556 , 0.00000393 ] )
      trcv_vs_za_fit_params[140] = numpy.array( [ 42.35740000 , -0.12224778 , 0.00070570 ] ) 
      trcv_vs_za_fit_params[160] = numpy.array( [ 42.44980000 , -0.06561444 , 0.00028481 ] ) 
      trcv_vs_za_fit_params[210] = numpy.array( [ 43.53470000 , -0.06032556 , 0.00094985 ] )
      trcv_vs_za_fit_params[220] = numpy.array( [ 46.33950000 , 0.07806778 , -0.00182215 ] )
      trcv_vs_za_fit_params[230] = numpy.array( [ 44.78490000 , 0.10587222 , -0.00187741 ] )
      trcv_vs_za_fit_params[280] = numpy.array( [ 50.56900000 , 0.08150778 , -0.00180215 ] )
      trcv_vs_za_fit_params[340] = numpy.array( [ 49.09320000 , -0.02795556 , 0.00065852 ] )
      trcv_vs_za_fit_params[345] = numpy.array( [ 51.47230000 , 0.01942444 , -0.00038193 ] )
      trcv_vs_za_fit_params[350] = numpy.array( [ 52.30620000 , 0.03021889 , -0.00063785 ] )
     
      if freq_mhz >= ( min(trcv_vs_za_fit_params.keys())-10 ) and freq_mhz <= ( max(trcv_vs_za_fit_params.keys())+10 ) :
         best_freq=None
         mindist=1e6
         for freq_key in trcv_vs_za_fit_params.keys() :
            dist = math.fabs( float(freq_key) - freq_mhz ) 
            if dist < mindist :
               mindist = dist
               best_freq = freq_key
         
         if best_freq is not None :
            fit_params = trcv_vs_za_fit_params[best_freq]            
            trcv_vs_za = fit_params[0] + fit_params[1]*za_deg + fit_params[2]*(za_deg**2)
            
            print("INFO : calculated trcv_aavs(%.2f MHz at za = %.2f [deg]) = %.2f [K] (vs. normal = %.2f [K]) , used fit coefficients %.4f, %.4f , %.4f" %  (freq_mhz,za_deg,trcv_vs_za,trcv,fit_params[0],fit_params[1],fit_params[2]))
            trcv = trcv_vs_za                           
         else :
            print("WARNING : best_freq not found for %.2f MHz -> do not using ZA dependence" % (freq_mhz))
      else :
         print("WARNING : frequency %.2f MHz outside the range of known ZA dependence (80-350 MHz) -> do not using ZA dependence" % (freq_mhz))
   else:
      print("WARNING : za_deg is None - using T_rcv(freq) only not as a function of pointing ZA[deg]")

   print("INFO : using trcv_aavs2( %.2f MHz) = %.2f [K] (cubic = %s)" % (freq_mhz,trcv,using_cubic))
      
   
   return trcv

# EDA2 - FEM contribution :
# see /home/msok/Desktop/EDA2/logbook/20200128_EDA2_receiver_temperature.odt
def trcv_eda2_fem( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too 
   print("INFO : using trcv_eda2 (same as EDA-1 + receiver temperature of the FEM in the SmartBox")

   x=[ 40, 50, 60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350 ]
   y=[ 13517.9, 13517.9,11947.6,11202.1,10830.2,10539,10205.4,10015.1,9683.08,9575.04,9396.75,9193.21,9046.74,8947.86,8820.11,8807.01,8596.01,8451,8367.9,8316.72,8230.98,8134.13,8023.28,7987.7,7761.26,7889.29,7796.44,7724.62,7741.32,7572.9,7504.12,7483.41 ]

#   gain_eda2_lna = 100.00 # about 20dB 
   gain_eda2_lna_db = eda2_lna_gain_db( freq_mhz )
   gain_eda2_lna = db2num( gain_eda2_lna_db )
   trcv_fem = 10000.00 / gain_eda2_lna # Approximately 10000 K / 20 dB of MWA LNA
   if freq_mhz >= 40 and freq_mhz <= 350 :
       trcv_fem_interpol = interp1d(x, y, kind='cubic')
       trcv_fem = float( trcv_fem_interpol( freq_mhz ) ) / gain_eda2_lna

   print("INFO : using trcv_eda2_fem( %.2f MHz) =  %.2f [K]" % (freq_mhz,trcv_fem))
   
   return trcv_fem

def interpol( x, x1, y1, x2, y2 ) :
   delta = 0.00
   if math.fabs( x2 - x1 ) > 0.000001  : 
      delta = ((y2-y1)/(x2-x1))*(x-x1)
   
   interpol_val = y1 + delta
   
   return interpol_val;


def interpolate( x_values, y_values, x ) :
   l = len(x_values)

   if x <= x_values[0] :
      return y_values[0]

   if x >= x_values[l-1] :
      return y_values[l-1]
      
   prev_idx = 0
   
   for i in range(1,l) :
      if math.fabs( x - x_values[i] ) < 0.00001 :
         return y_values[i]

      if x > x_values[prev_idx] and x < x_values[i] :
         return interpol( x, x_values[prev_idx], y_values[prev_idx], x_values[i], y_values[i] )
       
      prev_idx = i
   
   
   # error - not found
   return -1000

#
# sources :
#
# 50 MHz  :
# 70 MHz  : 
# 110 MHz :
# 160 MHz :
# 230 MHz : /home/msok/Desktop/SKA/papers/2021/EDA2_paper_Randall/20210602_Trcv_230MHz_data.odt, Fitted to all and mean : T_rcv_x ~ 133 K , T_rcv_y ~ 140 K 
# 320 MHz :
#
# see : /home/msok/Desktop/EDA2/papers/2021/EDA2_paper_Randall/20210518_trcv_recap_and_gain_vs_time.odt 
# /home/msok/Desktop/SKA/papers/2020/M_Sokolowski/Sensitivity_DB/20200903_EDA2_trcv_fitted_from_drift_scan_data.odt 
# /home/msok/Desktop/SKA/papers/2021/EDA2_paper_Randall/references/OLD_eda2_analysis/202009_update.odp
def trcv_eda2_single_dipole_fit( freq_mhz , use_cubic=False, pol="X", exact_dist=1.00, use_median_fit=False, use_eda2_paper=False ) :
   # FITTED AS IN : /home/msok/Desktop/SKA/papers/2020/M_Sokolowski/Sensitivity_DB/20200903_EDA2_trcv_fitted_from_drift_scan_data.odt
   # /home/msok/Desktop/EDA2/data/Trcv_fitted_from_the_sky/lst13-17hours
   # awk '{printf("%s , ",$1);}' eda2_fitted_trcv_lst13-17h_X.txt 
   # awk '{printf("%s , ",$3;;}' eda2_fitted_trcv_lst13-17h_X.txt 
   x_x_pol = [ 50 , 70 , 110 , 160 , 230 , 320 ]
   y_x_pol = [ 82760 , 3600 , 820 , 115 , 133 , 252 ]
 
   # /home/msok/Desktop/EDA2/data/Trcv_fitted_from_the_sky/lst13-17hours
   # awk '{printf("%s , ",$1);}' eda2_fitted_trcv_lst13-17h_Y.txt 
   # awk '{printf("%s , ",$3;;}' eda2_fitted_trcv_lst13-17h_Y.txt 
   x_y_pol = [ 50 , 70 , 110 , 160 , 230 , 320 ] 
   y_y_pol = [ 84630 , 2549 , 810 , 224 , 150 , 465 ]  # was [ 84630 , 2549 , 810 , 224 , 184 , 465 ] before : 20210602_Trcv_230MHz_data.odt etc 
   
   # fitted to median signle dipole lightcurves as in : /home/msok/Desktop/SKA/papers/2021/EDA2_paper_Randall/20210526_trcv_from_median_single_antenna_lc.odt  
   if use_median_fit :
      x_x_pol = [    50 ,   70 , 110 , 160 , 230 , 320 ]
      y_x_pol = [ 82760 , 3418 , 726 , 103 , 113 , 174 ] # TODO : fix 50 and 320 MHz value 
      
      x_y_pol = [    50 ,   70 , 110 , 160 , 230 , 320 ] 
      y_y_pol = [ 84630 , 2523 , 650 ,  79 ,  73 , 120 ] # TODO : fix 50 and 320 MHz value

   # As in EDA2 paper 2021 : /home/msok/Desktop/SKA/papers/2021/EDA2_paper_Randall/
   # 50 MHz  :
   # 70 MHz  : 
   # 110 MHz : 20200421_ch141_110MHz_analysis_FINAL.odt , 20200421_ch141_110MHz_analysis.odt , in 20210601_Trcv_110MHz_data.odt : T_rcv_x ~ 850 K , T_rcv_y ~ 850 K 
   # 160 MHz : 20200407_ch204_160MHz_analysis_FINAL.odt 20210603_Trcv_160MHz_data.odt : best fitting sensitivity are :  T_rcv_x=200 K and T_rcv_y = 160 K, but fitted : T_rcv_x ~ 110 K and T_rcv_y ~ 90 K 
   # 230 MHz : 20210602_Trcv_230MHz_data.odt, Fitted to all and mean : T_rcv_x ~ 133 K , T_rcv_y ~ 140 K 
   # 320 MHz : 20210601_Trcv_320MHz_data.odt Looks like T_rcv_x ~ 200 K , T_rcv_y ~ 250 K 
   if use_eda2_paper : 
      # updated according to /home/msok/Desktop/EDA2/papers/2021/EDA2_paper_Randall/20210601_paper_plots_INIT.odt and references 
      x_x_pol = [    50 ,   70 , 110 , 160 , 230 , 320 ]
      y_x_pol = [ 82760 , 1700 , 850 , 200 , 133 , 170 ] # was : [ 82760 , 3418 , 726 , 103 , 113 , 174 ]

      x_y_pol = [    50 ,   70 , 110 , 160 , 230 , 320 ] 
      y_y_pol = [ 84630 , 5000 , 850 , 160 , 150 , 235 ] # was : [ 84630 , 2523 , 650 , 160 , 150 , 235 ]


 
   l = len(x_x_pol)
 
   if pol == "X" :
      # try exact first +/- exact_dist to try to use the exact fitted values (not the interpolation results which can be strange ...)
      f_idx = 0
      for f_mhz in x_x_pol :
         if math.fabs( f_mhz - freq_mhz ) <= exact_dist :
            return y_x_pol[f_idx]
         
         f_idx += 1
   
      if freq_mhz >= x_x_pol[0] and freq_mhz <= x_x_pol[l-1] :
#         trcv_interpol = interp1d( x_x_pol, y_x_pol , kind="linear" ) # , kind='cubic')      
#         trcv_eda = trcv_interpol( freq_mhz )         
         trcv_eda = numpy.interp( freq_mhz, x_x_pol, y_x_pol, left=y_x_pol[0], right=y_x_pol[l-1] )
#         trcv_eda = interpolate( x_x_pol, y_x_pol, freq_mhz )
       
         return trcv_eda
      else :
         if freq_mhz <  x_x_pol[0] :
            return y_x_pol[0]
         else :
            return y_x_pol[l-1]
   elif pol == "Y" :
      # try exact first +/- exact_dist to try to use the exact fitted values (not the interpolation results which can be strange ...)
      f_idx = 0
      for f_mhz in x_y_pol :
         if math.fabs( f_mhz - freq_mhz ) <= exact_dist :
            return y_y_pol[f_idx]
         
         f_idx += 1


      if freq_mhz >= x_y_pol[0] and freq_mhz <= x_y_pol[l-1] :
#         trcv_interpol = interp1d( x_y_pol, y_y_pol , kind="linear" ) # , kind='cubic')
#         trcv_eda = trcv_interpol( freq_mhz )
         trcv_eda = numpy.interp( freq_mhz, x_y_pol, y_y_pol, left=y_y_pol[0], right=y_y_pol[l-1] )
#         trcv_eda = interpolate( x_y_pol, y_y_pol, freq_mhz )
       
         return trcv_eda
      else :
         if freq_mhz <  x_y_pol[0] :
            return y_y_pol[0]
         else :
            return y_y_pol[l-1]
          
   else :
      print("ERROR : polarisation %s is unknown" % (pol))
      
   # error :      
   return -1000   
    
 
 
   
# This is based on Daniel's paper : Figure 9 in https://arxiv.org/pdf/2003.05116.pdf
# Digitised : /home/msok/Desktop/MWA/papers/2020/Daniel/images/FINAL/Figure9_Trcv_vs_Freq_EDA2.png   
#             /home/msok/Desktop/MWA/papers/2020/Daniel/images/FINAL/Figure9_Trcv_vs_Freq_EDA2_SINGLE_ISOLATED_ELEMENT.txt
#             using plot_digitiser! 
def trcv_eda2_single_dipole( freq_mhz , use_cubic=False, add_fem=True ):   
   x = [ 55.6462, 57.2146, 57.8419, 58.7829, 60.0376, 60.3513, 61.6060, 62.5471, 63.4881, 64.4291, 65.0565, 66.6248, 67.8795, 68.5069, 69.1343, 70.3890, 71.6437, 72.4801, 73.5257, 74.5713, 75.6169, 76.4534, 77.4990, 79.3810, 80.4266, 81.6813, 82.7269, 83.3542, 84.1907, 85.4454, 86.7001, 87.9548, 89.4187, 91.0916, 92.7645, 94.0192, 95.2739, 97.1560, 98.6198, 100.084, 101.129, 102.593, 104.266, 105.939, 107.821, 109.076, 110.121, 111.167, 112.212, 113.363, 114.147, 115.245, 116.656, 118.381, 119.479, 120.577, 121.832, 123.087, 124.498, 126.380, 128.890, 129.987, 131.242, 132.497, 133.908, 135.320, 136.888, 138.143, 139.555, 140.809, 141.907, 143.162, 144.417, 145.985, 147.710, 148.965, 150.063, 151.317, 152.415, 153.513, 154.454, 155.709, 157.434, 159.473, 161.198, 162.923, 164.492, 167.001, 169.040, 170.295, 172.020, 173.118, 174.373, 175.627, 177.666, 179.705, 182.371, 184.253, 185.351, 187.233, 189.900, 190.997, 192.252, 193.036, 194.448, 195.546, 196.957, 197.898, 201.035, 203.492, 205.270, 206.420, 208.093, 210.079, 212.484, 214.262, 216.458, 217.817, 219.281, 220.117, 221.476, 223.777, 225.450, 226.495, 227.854, 228.795, 229.737, 231.200, 233.501, 235.383, 238.310, 241.133, 243.538, 245.002, 247.721, 250.230, 253.367, 255.144, 257.967, 260.477, 262.986, 265.182, 267.273, 270.201, 273.233, 275.115, 277.520, 280.761, 283.584, 287.453, 290.694, 293.726, 296.445, 298.223, 300.105, 302.300, 304.287, 305.646, 307.424, 309.515, 312.129, 314.638, 317.148, 319.448, 322.585, 324.780, 326.453 ]
   y = [ 9302.72, 7336.19, 6348.78, 5608.92, 4804.13, 4244.27, 3749.65, 3312.68, 2896.57, 2612.40, 2307.96, 2081.54, 1997.32, 1916.51, 1764.56, 1558.92, 1349.10, 1220.95, 1139.73, 1035.02, 939.927, 877.402, 813.418, 743.789, 694.312, 639.264, 592.646, 564.766, 549.427, 523.581, 505.866, 468.976, 437.780, 414.323, 394.832, 376.258, 361.034, 341.690, 325.616, 310.298, 295.700, 283.736, 266.691, 252.402, 240.528, 232.390, 224.528, 218.430, 209.592, 199.732, 191.651, 185.805, 177.369, 171.959, 167.577, 159.969, 155.090, 149.585, 142.794, 134.910, 126.806, 123.574, 118.574, 113.777, 109.173, 104.756, 100.518, 95.9539, 91.5975, 88.8036, 87.4388, 83.0391, 78.8608, 74.8927, 71.8625, 69.3118, 66.1649, 63.8165, 62.8357, 62.1902, 61.2344, 59.9828, 57.2595, 55.5130, 52.9926, 51.3762, 50.0670, 47.5478, 45.1553, 43.7780, 43.5525, 43.3283, 41.5752, 41.1481, 40.5157, 40.5157, 40.9362, 41.5752, 41.7904, 41.5752, 41.5752, 41.5752, 41.1481, 40.5157, 40.9362, 41.1481, 41.3611, 41.5752, 41.5752, 42.0067, 41.7185, 41.7185, 43.0310, 44.3848, 44.2323, 45.0001, 45.0001, 45.4672, 46.2564, 47.2216, 47.5478, 47.5478, 47.3844, 47.2216, 47.2216, 47.0593, 47.5478, 48.5399, 49.8950, 50.2396, 49.7235, 49.8950, 50.0670, 52.3581, 54.5658, 55.1322, 54.5658, 54.3783, 54.9427, 56.2825, 56.8667, 57.6550, 58.6558, 59.2646, 58.4543, 58.4543, 59.0610, 59.6739, 60.9191, 60.5012, 61.3399, 62.4046, 62.4046, 62.4046, 62.4046, 61.5514, 61.5514, 61.5514, 60.9191, 60.0861, 60.5012, 60.2933, 60.0861, 61.3399, 61.9765, 63.0523, 63.9264 ]
   
   trcv_eda_lna = 9302.72
   if freq_mhz >= 326 :
      trcv_eda_lna = 63.9264
      
   if freq_mhz >= x[0] and freq_mhz <= x[len(x)-1] :
      trcv_interpol = interp1d(x, y, kind='cubic')
      trcv_eda_lna = float( trcv_interpol( freq_mhz ) )
      
   trcv_out = trcv_eda_lna
   
   trcv_fem = 0.00
   if add_fem :   
      trcv_fem = trcv_eda2_fem( freq_mhz )
      print("DEBIG : %.4f (%.f) + %.4f = ???" % (trcv_out,trcv_eda_lna,trcv_fem))     
      trcv_out += trcv_fem

   print("INFO : using trcv_eda2_single_dipole( %.2f MHz) = %.2f + %.2f  = %.2f [K]" % (freq_mhz,trcv_eda_lna,trcv_fem,trcv_out))
      
   return trcv_out

# This is based on Daniel's paper : Figure 9 in https://arxiv.org/pdf/2003.05116.pdf
# Digitised : /home/msok/Desktop/MWA/papers/2020/Daniel/images/FINAL/Figure9_Trcv_vs_Freq_EDA2.png   
#             /home/msok/Desktop/MWA/papers/2020/Daniel/images/FINAL/Figure9_Trcv_vs_Freq_EDA2_SINGLE_ISOLATED_ELEMENT.txt
#             using plot_digitiser! 
def trcv_eda1_single_dipole( freq_mhz , use_cubic=False, add_fem=True ):
   ret = trcv_eda2_single_dipole( freq_mhz , use_cubic=use_cubic, add_fem=False ) # add_fem=False -> EDA1 Trcv as in the Daniel's paper for a single EDA2 dipole 
   
   return ret   



# EDA2 : 
# see /home/msok/Desktop/EDA2/logbook/20200128_EDA2_receiver_temperature.odt
def trcv_eda2( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too 
   print("INFO : using trcv_eda2 (same as EDA-1 + receiver temperature of the FEM in the SmartBox")

   x=[      40,      50,      60,      70,      80,    90,     100,     110, 120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350 ]
   y=[ 13517.9, 13517.9, 11947.6, 11202.1, 10830.2, 10539, 10205.4, 10015.1, 9683.08,9575.04,9396.75,9193.21,9046.74,8947.86,8820.11,8807.01,8596.01,8451,8367.9,8316.72,8230.98,8134.13,8023.28,7987.7,7761.26,7889.29,7796.44,7724.62,7741.32,7572.9,7504.12,7483.41 ]

#   gain_eda2_lna = 100.00 # about 20dB 
   gain_eda2_lna_db = eda2_lna_gain_db( freq_mhz )
   gain_eda2_lna = db2num( gain_eda2_lna_db )
   trcv_fem = 10000.00 / gain_eda2_lna # Approximately 10000 K / 20 dB of MWA LNA
   if freq_mhz >= 40 and freq_mhz <= 350 :
       trcv_fem_interpol = interp1d(x, y, kind='cubic')
       trcv_fem = float( trcv_fem_interpol( freq_mhz ) ) / gain_eda2_lna

   # at the moment return the same as for EDA-1 (or EDA-2)
   trcv_eda_lna = trcv_from_skymodel_with_err( freq_mhz , use_cubic )
   trcv_out = trcv_eda_lna + trcv_fem
   print("INFO : using trcv_eda2( %.2f MHz) = %.2f + %.2f  = %.2f [K]" % (freq_mhz,trcv_eda_lna,trcv_fem,trcv_out))
   
   return trcv_out


# EDA1 (as for the EDA1 paper) :
def trcv_eda1( freq_mhz , use_cubic=False ): # default was power law fit but I wanted to be able to use cubic fit too
   trcv_out = trcv_from_skymodel_with_err( freq_mhz , use_cubic )
   
   print("INFO : using trcv_eda1( %.2f MHz) = %.2f [K]" % (freq_mhz,trcv_out))
   
   return trcv_out


def trcv_from_skymodel_with_err_cubic( freq_mhz ) :
    x=[55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0]
    y=[3722.55055596,1421.83480972,655.8033454,409.69506365,356.59552554,293.4326318,200.26593869,163.91940732,157.46121101,112.19530177,51.14849704,68.56693587,28.33318976,43.72596282,34.41633474,62.76422256]
    tlna_cubic = interp1d(x, y, kind='cubic')
    trcv = float( tlna_cubic(freq_mhz) )
    
    return trcv


def trcv_lightcurve201612( freq_mhz ):
    T_rcv =  +53480.71172570309136062860488891601562500000000000000000 - 1455.76561026112926811038050800561904907226562500000000*freq_mhz \
             +12.16264789947674174186431628186255693435668945312500*freq_mhz*freq_mhz -0.00271098517125182300949171043669139180565252900124*freq_mhz*freq_mhz*freq_mhz \
             -0.00030428122418899071835451941581140999915078282356*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000046193435183846278872557794374642536894270961*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000959650510252193562322564293230395260358989*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000001953073405836920809493335023823663365930*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             -0.00000000000028926757335617184960944065203237193611*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
             +0.00000000000000055485231182041931661150443032917814*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz
             
    return T_rcv             
   
def trcv_fit_lightcurve_20160703_polfit( freq_mhz ):
   t_lna =  +36853.38919790466752601787447929382324218750000000000000  -1377.64458069924739902489818632602691650390625000000000 * freq_mhz +20.37196908756931890138730523176491260528564453125000 * freq_mhz*freq_mhz \
            -0.13515251806813355361924777753301896154880523681641 * freq_mhz*freq_mhz*freq_mhz +0.00022215700342832035230414278181854115246096625924 * freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            +0.00000169842384409759492976880849990362065682347747 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz -0.00000000464297999395697167928294570157440879221156 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            -0.00000000003336389317067374617426693716906293402102 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz +0.00000000000017598078793649799250229898970650520328 * freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz \
            -0.00000000000000023477806740472191603243644100388291 * 1*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz*freq_mhz

   return t_lna


def trcv_multi( freq_mhz , type, use_cubic=False, za_deg=0, pol="X", use_median_fit=False, use_eda2_paper=False ):
   t_rcv = 0.00
   if type.lower() == "eda2" or type.lower() == "trcv_eda2" :
      t_rcv = trcv_eda2( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "eda1" or type.lower() == "trcv_eda1" :
      t_rcv = trcv_eda1( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "aavs2" or type.lower() == "trcv_aavs2" :
      t_rcv = trcv_aavs2( freq_mhz, use_cubic=use_cubic )       
   elif type.lower() == "eda2_singledip" or type.lower() == "trcv_eda2_singledip" or type.lower() == "trcv_eda2_single_dipole" or type.lower() == "trcv_eda2_single_dip" :
      t_rcv = trcv_eda2_single_dipole( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "eda1_singledip" or type.lower() == "trcv_eda1_singledip" or type.lower() == "trcv_eda1_single_dipole" or type.lower() == "trcv_eda1_single_dip" :
      t_rcv = trcv_eda1_single_dipole( freq_mhz, use_cubic=use_cubic )
   elif type.lower() == "trcv_aavs2_vs_za_deg" :
      t_rcv = trcv_aavs2( freq_mhz, za_deg=za_deg )
   elif type.lower() == "trcv_eda2_single_dipole_fit" :      
      if pol == "X" or pol == "x" :
         t_rcv = trcv_eda2_single_dipole_fit( freq_mhz, use_cubic=use_cubic, pol="X", use_median_fit=use_median_fit, use_eda2_paper=use_eda2_paper )
      elif pol == "Y" or pol == "y" :
         t_rcv = trcv_eda2_single_dipole_fit( freq_mhz, use_cubic=use_cubic, pol="Y", use_median_fit=use_median_fit, use_eda2_paper=use_eda2_paper )
      else :
         print("ERROR : unknown polarisation %s -> crashin script now (error in code)" % (pol))
         sys,exit(-1)
   elif type.lower() == "trcv_eda2_single_dipole_fit_x" :
      t_rcv = trcv_eda2_single_dipole_fit( freq_mhz, use_cubic=use_cubic, pol="X", use_median_fit=use_median_fit, use_eda2_paper=use_eda2_paper )
   elif type.lower() == "trcv_eda2_single_dipole_fit_y" :
      t_rcv = trcv_eda2_single_dipole_fit( freq_mhz, use_cubic=use_cubic, pol="Y", use_median_fit=use_median_fit, use_eda2_paper=use_eda2_paper )
   else :
      print("ERROR : unknown receiver temperature type = %s -> CRITICAL ERROR -> cannot continue" % (type))
      exit(-1)
      # t_rcv = trcv_from_skymodel_with_err_cubic( freq_mhz )

   return t_rcv


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tCalculate receiver temperature (T_rcv) for different SKA-Low stations and configurations\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-t','--type','--trcv_type',dest="trcv_type",default="eda2",help="Type of receiver temperature (station/configuration) [default %default]",metavar="STRING")
   parser.add_option('-s','--start_freq','--start_freq_mhz',dest="start_freq_mhz",default=160,help="Start of frequency range to calculate T_rcv [default %default]",type="float")
   parser.add_option('-e','--end_freq','--end_freq_mhz',dest="end_freq_mhz",default=160,help="End of frequency range to calculate T_rcv [default %default]",type="float")
   parser.add_option('-m','--use_median_fit','--trcv_fitted_to_median_curve','--trcv_from_median_curve','--trcv_from_median_lc',action="store_true",dest="use_median_fit",default=False, help="Use receiver temperature fitted to median lightcurve of all antennas [default %default]")
   parser.add_option('-p','--eda2_paper','--use_eda2_paper','--eda2_2021_paper','--paper',action="store_true",dest="use_eda2_paper",default=False, help="Use receiver temperature exactly as used in the calculation of station sensitivity for the 2021 EDA2 paper [default %default]")

   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)


if __name__ == "__main__":
    (options, args) = parse_options()

#    start_freq_mhz = 160.00
#    if len(sys.argv) >= 1:
#       start_freq_mhz=float( sys.argv[1] )

#    end_freq_mhz = 160.00
#    if len(sys.argv) >= 2:
#       end_freq_mhz=float( sys.argv[2] )


#    type="eda2"
#    if len(sys.argv) >= 3:
#       type=sys.argv[3]

    step_mhz = 1.00
    freq_mhz = options.start_freq_mhz
    while freq_mhz <= options.end_freq_mhz :
       t_rcv = trcv_multi( freq_mhz , options.trcv_type, False, use_median_fit=options.use_median_fit, use_eda2_paper=options.use_eda2_paper )
       print("%s : T_rcv ( %.2f MHz) = %.8f [K]" % (options.trcv_type,freq_mhz,t_rcv))
       
       freq_mhz += step_mhz
       
       
    
    