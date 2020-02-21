# Mapping of TPMs to antennas in AAVS :

class AntennaLocation :
   x    = -1e6
   y    = -1e6
   base = -1e6
   tpm  = -1   
   
   def __init__(self, x, y, base, tpm):
        self.x = x
        self.y = y
        self.base = base
        self.tpm  = tpm
        
   def __repr__(self):
        out_line = "(%.2f,%.2f) base = %d tpm = %d" % (self.x,self.y,self.base,self.tpm)
        return out_line 
        
   def __str__(self):
        out_line = "Antenna at (x,y) = (%.2f,%.2f) , base = %d, tpm = %d" % (self.x,self.y,self.base,self.tpm)
        return out_line 


class Station :
    tpm_list = {}

    def __init__(self) :
        self.init_tpm_list()

    
    def init_tpm_list( self ) :
        self.tpm_list = {}

        self.tpm_list['1'] = [ AntennaLocation( 8.1750, 1.0600, 50, 1), AntennaLocation( 5.3560, 1.0780, 83, 1), AntennaLocation( 10.1610, 1.9410, 33, 1), AntennaLocation( 10.2520, 3.5290, 32, 1), AntennaLocation( 13.7700, 3.2680, 22, 1), AntennaLocation( 15.7210, 4.8030, 8, 1), AntennaLocation( 12.1720, 2.2730, 21, 1), AntennaLocation( 6.7930, 1.2800, 64, 1), AntennaLocation( 2.4050, -5.8730, 103, 1), AntennaLocation( 2.5310, 1.6440, 99, 1), AntennaLocation( 0.8650, -2.5880, 119, 1), AntennaLocation( 2.4170, -2.4290, 101, 1), AntennaLocation( 7.7290, -11.3110, 70, 1), AntennaLocation( 7.4410, -14.7800, 72, 1), AntennaLocation( 2.6080, -7.4060, 104, 1),  ]
        self.tpm_list['2'] = [ AntennaLocation( 4.3460, 3.3240, 84, 2), AntennaLocation( 14.7050, 6.1960, 7, 2), AntennaLocation( 13.7410, 4.8470, 23, 2), AntennaLocation( 11.3590, 8.9840, 29, 2), AntennaLocation( 14.6100, 8.2100, 6, 2), AntennaLocation( 13.8230, 10.8530, 25, 2), AntennaLocation( 12.8060, 8.2660, 24, 2), AntennaLocation( 11.2500, 6.6880, 30, 2), AntennaLocation( 9.1090, 5.3370, 51, 2), AntennaLocation( 9.8970, 7.1540, 53, 2), AntennaLocation( 12.0730, 11.9410, 26, 2), AntennaLocation( 8.2760, 6.5910, 52, 2), AntennaLocation( 11.2980, 4.5760, 31, 2), AntennaLocation( 7.0000, 3.0000, 63, 2), AntennaLocation( 6.4320, 4.9460, 61, 2), AntennaLocation( 7.8090, 4.8870, 62, 2),  ]
        self.tpm_list['3'] = [ AntennaLocation( 11.5020, 13.2030, 27, 3), AntennaLocation( 6.1760, 9.9870, 59, 3), AntennaLocation( 6.9680, 16.0510, 58, 3), AntennaLocation( 8.6620, 12.2690, 56, 3), AntennaLocation( 8.5700, 13.8310, 57, 3), AntennaLocation( 4.2300, 5.9670, 85, 3), AntennaLocation( 9.9940, 10.6600, 55, 3), AntennaLocation( 5.6460, 6.3950, 86, 3), AntennaLocation( 10.4600, 12.2070, 28, 3), AntennaLocation( 3.6440, 4.6650, 98, 3), AntennaLocation( 4.7830, 10.9620, 88, 3), AntennaLocation( 5.0830, 8.3120, 87, 3), AntennaLocation( 5.9970, 12.1520, 89, 3), AntennaLocation( 6.7680, 7.4230, 60, 3), AntennaLocation( 5.7260, 14.6650, 91, 3), AntennaLocation( 8.8240, 9.2730, 54, 3),  ]
        self.tpm_list['4'] = [ AntennaLocation( 4.3450, 17.0720, 92, 4), AntennaLocation( 4.1010, 14.3110, 90, 4), AntennaLocation( 2.9220, 12.9280, 94, 4), AntennaLocation( 0.5890, 9.5060, 125, 4), AntennaLocation( 0.9060, 6.0910, 123, 4), AntennaLocation( 1.4250, 11.8290, 126, 4), AntennaLocation( 0.4430, 1.7010, 120, 4), AntennaLocation( 1.8660, 8.9920, 124, 4), AntennaLocation( 1.0780, 2.9630, 121, 4), AntennaLocation( 3.0840, 10.1830, 95, 4), AntennaLocation( 1.6440, 16.9370, 129, 4), AntennaLocation( 2.4960, 15.3440, 93, 4), AntennaLocation( 2.8380, 6.3840, 97, 4), AntennaLocation( 3.3020, 8.6000, 96, 4), AntennaLocation( 1.3050, 13.6070, 128, 4), AntennaLocation( 1.7270, 4.2670, 122, 4),  ]
        self.tpm_list['5'] = [ AntennaLocation( -2.6420, 8.3150, 156, 5), AntennaLocation( -0.7600, 13.7450, 131, 5), AntennaLocation( -2.1520, 13.1770, 159, 5), AntennaLocation( -2.2170, 15.3140, 161, 5), AntennaLocation( -2.8660, 17.3950, 162, 5), AntennaLocation( -3.6660, 14.2790, 160, 5), AntennaLocation( -2.0490, 10.6690, 157, 5), AntennaLocation( 0.0010, 12.5450, 127, 5), AntennaLocation( -4.7480, 16.0280, 163, 5), AntennaLocation( -1.0790, 6.4650, 135, 5), AntennaLocation( -0.3590, 3.2620, 136, 5), AntennaLocation( -3.2150, 12.0950, 158, 5), AntennaLocation( -0.1440, 10.9170, 132, 5), AntennaLocation( -0.5120, 8.2030, 134, 5), AntennaLocation( -0.3160, 16.1120, 130, 5), AntennaLocation( -0.8420, 9.5510, 133, 5),  ]
        self.tpm_list['6'] = [ AntennaLocation( -5.3870, 10.7130, 165, 6), AntennaLocation( -9.2320, 10.2790, 195, 6), AntennaLocation( -5.3130, 8.1110, 167, 6), AntennaLocation( -7.5630, 10.3890, 190, 6), AntennaLocation( -8.3460, 11.7230, 194, 6), AntennaLocation( -6.8060, 13.8320, 191, 6), AntennaLocation( -2.8330, 3.6360, 153, 6), AntennaLocation( -4.4820, 9.3780, 166, 6), AntennaLocation( -6.4010, 6.7790, 188, 6), AntennaLocation( -2.1250, 5.5330, 155, 6), AntennaLocation( -3.9900, 5.5160, 154, 6), AntennaLocation( -9.7510, 12.9710, 193, 6), AntennaLocation( -12.2530, 12.3530, 224, 6), AntennaLocation( -7.8800, 15.0190, 192, 6), AntennaLocation( -5.4940, 12.2140, 164, 6), AntennaLocation( -6.5710, 9.3240, 189, 6),  ]
        self.tpm_list['7'] = [ AntennaLocation( -12.5090, 9.7610, 226, 7), AntennaLocation( -13.3740, 11.3320, 225, 7), AntennaLocation( -8.1760, 8.1180, 197, 7), AntennaLocation( -12.1020, 8.1090, 227, 7), AntennaLocation( -4.3900, 2.6470, 168, 7), AntennaLocation( -16.0750, 6.9830, 251, 7), AntennaLocation( -9.4280, 7.3090, 198, 7), AntennaLocation( -10.8460, 6.0360, 221, 7), AntennaLocation( -11.2360, 11.3250, 223, 7), AntennaLocation( -14.8820, 9.4440, 250, 7), AntennaLocation( -9.2520, 4.4680, 200, 7), AntennaLocation( -10.8450, 9.9520, 222, 7), AntennaLocation( -8.3320, 6.2870, 199, 7), AntennaLocation( -7.0340, 5.4000, 187, 7), AntennaLocation( -9.9110, 8.8870, 196, 7), AntennaLocation( -14.2880, 6.4610, 249, 7),  ]
        self.tpm_list['8'] = [ AntennaLocation( -15.5190, 1.3100, 245, 8), AntennaLocation( -17.3270, 3.1210, 252, 8), AntennaLocation( -15.7080, 2.7420, 246, 8), AntennaLocation( -8.5860, 3.0280, 201, 8), AntennaLocation( -11.9210, 3.8230, 220, 8), AntennaLocation( -6.7210, 2.5290, 186, 8), AntennaLocation( -15.8550, 5.1710, 248, 8), AntennaLocation( -13.6270, 0.6210, 230, 8), AntennaLocation( -14.9250, 3.9810, 247, 8), AntennaLocation( -5.8390, 0.3280, 169, 8), AntennaLocation( -10.5430, 2.9360, 219, 8), AntennaLocation( -16.9920, 1.3100, 253, 8), AntennaLocation( -9.2730, 0.9120, 202, 8), AntennaLocation( -12.9530, 2.7800, 229, 8), AntennaLocation( -13.0490, 4.6150, 228, 8), AntennaLocation( -11.0140, 1.0580, 218, 8),  ]
        self.tpm_list['9'] = [ AntennaLocation( 15.0260, -1.8480, 10, 9), AntennaLocation( 14.7130, -3.5750, 11, 9), AntennaLocation( 16.4280, -3.6580, 2, 9), AntennaLocation( 11.0710, -3.6340, 37, 9), AntennaLocation( 14.0760, -0.6410, 9, 9), AntennaLocation( 12.6020, 0.1220, 19, 9), AntennaLocation( 16.6260, 3.0440, 5, 9), AntennaLocation( 12.6690, -1.8380, 18, 9), AntennaLocation( 11.1740, -2.1270, 36, 9), AntennaLocation( 6.6140, -2.7280, 66, 9), AntennaLocation( 8.1950, -4.6060, 47, 9), AntennaLocation( 12.5510, -7.2090, 16, 9), AntennaLocation( 16.2080, -6.3200, 1, 9), AntennaLocation( 16.8070, 0.9430, 4, 9), AntennaLocation( 10.6490, -0.7710, 35, 9), AntennaLocation( 11.5180, 0.9910, 34, 9),  ]
        self.tpm_list['10'] = [ AntennaLocation( 4.9330, -6.8040, 79, 10), AntennaLocation( 6.1220, -8.2430, 68, 10), AntennaLocation( 9.4380, -10.8390, 43, 10), AntennaLocation( 4.5150, -8.2340, 78, 10), AntennaLocation( 9.7090, -9.3860, 44, 10), AntennaLocation( 5.1640, -9.8770, 77, 10), AntennaLocation( 13.5780, -9.4580, 15, 10), AntennaLocation( 2.9470, -4.1170, 102, 10), AntennaLocation( 10.9720, -6.5110, 38, 10), AntennaLocation( 14.1710, -5.1100, 12, 10), AntennaLocation( 14.2170, -7.2680, 13, 10), AntennaLocation( 9.3340, -2.7290, 48, 10), AntennaLocation( 12.5600, -4.0880, 17, 10), AntennaLocation( 6.9350, -1.0440, 65, 10), AntennaLocation( 13.8770, 0.8310, 20, 10), AntennaLocation( 9.1510, -0.4540, 49, 10),  ]
        self.tpm_list['11'] = [ AntennaLocation( 11.1120, -9.6780, 39, 11), AntennaLocation( 5.2900, -5.3650, 80, 11), AntennaLocation( 9.9150, -7.7940, 46, 11), AntennaLocation( 12.3160, -11.0720, 14, 11), AntennaLocation( 7.8400, -9.8770, 69, 11), AntennaLocation( 8.4110, -8.3500, 45, 11), AntennaLocation( 4.9570, -1.8330, 81, 11), AntennaLocation( 7.3210, -5.7250, 67, 11), AntennaLocation( 2.1780, -0.6590, 100, 11), AntennaLocation( 4.9520, -11.4300, 76, 11), AntennaLocation( 5.6850, -15.6200, 74, 11), AntennaLocation( 6.8110, -12.7670, 71, 11), AntennaLocation( 8.4530, -12.5050, 42, 11), AntennaLocation( 10.1060, -12.4750, 40, 11), AntennaLocation( 4.3900, -0.1320, 82, 11),  ]
        self.tpm_list['12'] = [ AntennaLocation( 1.9020, -11.2430, 115, 12), AntennaLocation( 1.4210, -13.9180, 112, 12), AntennaLocation( 1.2630, -17.5540, 110, 12), AntennaLocation( 2.7640, -8.8560, 105, 12), AntennaLocation( 0.0560, -13.5070, 113, 12), AntennaLocation( 0.9870, -5.5020, 117, 12), AntennaLocation( 2.8980, -14.4180, 107, 12), AntennaLocation( 3.6270, -11.9370, 106, 12), AntennaLocation( 2.6710, -15.9200, 108, 12), AntennaLocation( 0.3040, -11.4400, 114, 12), AntennaLocation( 0.5070, -15.7350, 111, 12), AntennaLocation( 2.8320, -17.2910, 109, 12), AntennaLocation( 4.0720, -13.2560, 75, 12), AntennaLocation( 0.1210, -7.7790, 116, 12), AntennaLocation( 4.6390, -16.8680, 73, 12), AntennaLocation( 0.0840, -4.2700, 118, 12),  ]
        self.tpm_list['13'] = [ AntennaLocation( -2.3150, -14.0160, 143, 13), AntennaLocation( -1.4310, -10.0010, 140, 13), AntennaLocation( -4.4830, -9.0850, 173, 13), AntennaLocation( -2.7920, -15.6380, 142, 13), AntennaLocation( -4.1390, -13.9780, 176, 13), AntennaLocation( -1.5420, -4.7980, 138, 13), AntennaLocation( -4.8340, -11.4530, 174, 13), AntennaLocation( -5.4990, -12.8450, 175, 13), AntennaLocation( -4.1390, -16.3830, 177, 13), AntennaLocation( -3.0800, -11.0450, 145, 13), AntennaLocation( -3.4750, -12.6980, 144, 13), AntennaLocation( -2.1260, -7.7130, 147, 13), AntennaLocation( -2.4850, -9.0920, 146, 13), AntennaLocation( -0.9640, -6.5050, 139, 13), AntennaLocation( -6.4590, -15.2010, 178, 13), AntennaLocation( -1.4500, -17.0410, 141, 13),  ]
        self.tpm_list['14'] = [ AntennaLocation( -10.7880, -9.3240, 213, 14), AntennaLocation( -11.7120, -10.7410, 212, 14), AntennaLocation( -9.3110, -9.2780, 207, 14), AntennaLocation( -2.9050, -5.7850, 149, 14), AntennaLocation( -13.9220, -10.6490, 235, 14), AntennaLocation( -8.7780, -10.8150, 208, 14), AntennaLocation( -6.8140, -11.8300, 180, 14), AntennaLocation( -10.1750, -13.6690, 210, 14), AntennaLocation( -12.3320, -12.5100, 236, 14), AntennaLocation( -6.0050, -9.0210, 182, 14), AntennaLocation( -8.8930, -7.6510, 206, 14), AntennaLocation( -3.9300, -7.6760, 148, 14), AntennaLocation( -7.4820, -9.7460, 181, 14), AntennaLocation( -10.2440, -11.6760, 211, 14), AntennaLocation( -9.0240, -12.8890, 209, 14), AntennaLocation( -7.8040, -14.4660, 179, 14),  ]
        self.tpm_list['15'] = [ AntennaLocation( -14.2030, -8.9540, 237, 15), AntennaLocation( -13.0520, -4.3920, 232, 15), AntennaLocation( -16.3280, -5.9130, 256, 15), AntennaLocation( -5.7680, -3.3080, 171, 15), AntennaLocation( -6.8890, -5.1750, 183, 15), AntennaLocation( -8.2100, -4.0850, 204, 15), AntennaLocation( -14.7800, -7.6570, 238, 15), AntennaLocation( -1.0250, -1.5030, 137, 15), AntennaLocation( -12.8900, -7.2080, 233, 15), AntennaLocation( -2.4280, -3.1710, 150, 15), AntennaLocation( -8.3740, -5.6100, 205, 15), AntennaLocation( -10.6360, -6.9570, 214, 15), AntennaLocation( -14.1520, -5.8720, 239, 15), AntennaLocation( -11.8970, -5.4980, 215, 15), AntennaLocation( -4.7300, -4.8080, 172, 15), AntennaLocation( -12.5750, -8.7410, 234, 15),  ]
        self.tpm_list['16'] = [ AntennaLocation( -12.3770, -0.7440, 231, 16), AntennaLocation( -15.9000, -1.7680, 242, 16), AntennaLocation( -3.1950, -1.3640, 151, 16), AntennaLocation( -15.2570, -0.4190, 244, 16), AntennaLocation( -14.0790, -1.4170, 243, 16), AntennaLocation( -9.1660, -1.8240, 203, 16), AntennaLocation( -10.0460, -3.6410, 216, 16), AntennaLocation( -2.8060, 1.0110, 152, 16), AntennaLocation( -4.7760, -0.8150, 170, 16), AntennaLocation( -7.4370, -0.8020, 185, 16), AntennaLocation( -6.4170, -1.8630, 184, 16), AntennaLocation( -11.1860, -2.1740, 217, 16), AntennaLocation( -14.1980, -3.4040, 241, 16), AntennaLocation( -16.9430, -0.3060, 254, 16), AntennaLocation( -15.1650, -4.6200, 240, 16), AntennaLocation( -16.0990, -3.1340, 255, 16),  ]

        return self.tpm_list


    def get_xy_pos( self, tpm_request=None ) :
        xs = []
        ys = []
        
        for t in range(1,len(self.tpm_list)+1) :
            if tpm_request is None or t == tpm_request :
                tpm = self.tpm_list[str(t)]   
                
                for ant in tpm :            
                   xs.append( ant.x )
                   ys.append( ant.y )
                
        return (xs,ys)
        
                