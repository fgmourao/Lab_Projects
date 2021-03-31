EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title "HeadstageNNC.sch"
Date "27 FEV 2021"
Rev "1"
Comp "UFMG"
Comment1 "Núcleo de Neurociências"
Comment2 "Márcio Flávio Dutra Moraes"
Comment3 "Paulo Aparecido Amaral Júnior"
Comment4 "Flávio Afonso Gonçalves Mourão"
$EndDescr
$Comp
L RHD2132_headstage-rescue:C_0402 C1
U 1 1 57309091
P 1200 6200
F 0 "C1" H 1260 6214 70  0000 L BNN
F 1 "100n" H 1259 6015 70  0000 L BNN
F 2 "Resistors_SMD:R_0603" H 1190 5990 65  0001 L TNN
F 3 "" H 1200 6200 60  0001 C CNN
	1    1200 6200
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:C_0402 C2
U 1 1 57309159
P 5125 6075
F 0 "C2" H 5185 6089 70  0000 L BNN
F 1 "10n" H 5184 5890 70  0000 L BNN
F 2 "Resistors_SMD:R_0603" H 5115 5865 65  0001 L TNN
F 3 "" H 5125 6075 60  0001 C CNN
	1    5125 6075
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:LETTER_L #FRAME1
U 1 1 57309285
P 2025 7775
F 0 "#FRAME1" H 2075 7825 60  0001 C CNN
F 1 "~" H 2025 -14825 60  0001 C CNN
F 2 "" H 2025 -14825 60  0001 C CNN
F 3 "" H 2025 -14825 60  0001 C CNN
	1    2025 7775
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND01
U 1 1 57309415
P 2925 4875
F 0 "#GND01" H 2826 4775 70  0001 L BNN
F 1 "GND" V 2925 4650 70  0000 L BNN
F 2 "" H 2925 4875 60  0001 C CNN
F 3 "" H 2925 4875 60  0001 C CNN
	1    2925 4875
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND02
U 1 1 57309479
P 2925 4975
F 0 "#GND02" H 2826 4875 70  0001 L BNN
F 1 "GND" V 2950 4750 70  0000 L BNN
F 2 "" H 2925 4975 60  0001 C CNN
F 3 "" H 2925 4975 60  0001 C CNN
	1    2925 4975
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND03
U 1 1 573094DD
P 5825 5575
F 0 "#GND03" H 5725 5475 70  0001 L BNN
F 1 "GND" V 5850 5350 70  0000 L BNN
F 2 "" H 5825 5575 60  0001 C CNN
F 3 "" H 5825 5575 60  0001 C CNN
	1    5825 5575
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND04
U 1 1 57309541
P 5825 5375
F 0 "#GND04" H 5725 5275 70  0001 L BNN
F 1 "GND" V 5850 5150 70  0000 L BNN
F 2 "" H 5825 5375 60  0001 C CNN
F 3 "" H 5825 5375 60  0001 C CNN
	1    5825 5375
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND05
U 1 1 573095A5
P 3725 6075
F 0 "#GND05" H 3625 5975 70  0001 L BNN
F 1 "GND" H 3625 5975 70  0000 L BNN
F 2 "" H 3725 6075 60  0001 C CNN
F 3 "" H 3725 6075 60  0001 C CNN
	1    3725 6075
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND06
U 1 1 57309609
P 5825 4775
F 0 "#GND06" H 5725 4675 70  0001 L BNN
F 1 "GND" V 5850 4550 70  0000 L BNN
F 2 "" H 5825 4775 60  0001 C CNN
F 3 "" H 5825 4775 60  0001 C CNN
	1    5825 4775
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND07
U 1 1 57309735
P 5375 6875
F 0 "#GND07" H 5275 6775 70  0001 L BNN
F 1 "GND" V 5400 6650 70  0000 L BNN
F 2 "" H 5375 6875 60  0001 C CNN
F 3 "" H 5375 6875 60  0001 C CNN
	1    5375 6875
	-1   0    0    1   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND08
U 1 1 57309799
P 1200 6600
F 0 "#GND08" H 1100 6500 70  0001 L BNN
F 1 "GND" H 1100 6500 70  0000 L BNN
F 2 "" H 1200 6600 60  0001 C CNN
F 3 "" H 1200 6600 60  0001 C CNN
	1    1200 6600
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND09
U 1 1 573097FD
P 5125 6475
F 0 "#GND09" H 5025 6375 70  0001 L BNN
F 1 "GND" H 5025 6375 70  0000 L BNN
F 2 "" H 5125 6475 60  0001 C CNN
F 3 "" H 5125 6475 60  0001 C CNN
	1    5125 6475
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+010
U 1 1 57309C49
P 5825 4975
F 0 "#P+010" V 5726 4875 70  0001 L BNN
F 1 "VCC" V 5850 4975 70  0000 L BNN
F 2 "" H 5825 4975 60  0001 C CNN
F 3 "" H 5825 4975 60  0001 C CNN
	1    5825 4975
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+011
U 1 1 57309D11
P 4825 6075
F 0 "#P+011" V 4725 5975 70  0001 L BNN
F 1 "VCC" V 4864 6102 70  0000 L BNN
F 2 "" H 4825 6075 60  0001 C CNN
F 3 "" H 4825 6075 60  0001 C CNN
	1    4825 6075
	-1   0    0    1   
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+012
U 1 1 57309D75
P 2925 5075
F 0 "#P+012" V 2825 4975 70  0001 L BNN
F 1 "VCC" V 2950 5100 70  0000 L BNN
F 2 "" H 2925 5075 60  0001 C CNN
F 3 "" H 2925 5075 60  0001 C CNN
	1    2925 5075
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+013
U 1 1 57309DD9
P 1200 5900
F 0 "#P+013" V 1100 5800 70  0001 L BNN
F 1 "VCC" H 1275 5850 70  0000 L BNN
F 2 "" H 1200 5900 60  0001 C CNN
F 3 "" H 1200 5900 60  0001 C CNN
	1    1200 5900
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+014
U 1 1 57309E3D
P 5275 6875
F 0 "#P+014" V 5175 6775 70  0001 L BNN
F 1 "VCC" V 5350 6900 70  0000 L BNN
F 2 "" H 5275 6875 60  0001 C CNN
F 3 "" H 5275 6875 60  0001 C CNN
	1    5275 6875
	-1   0    0    -1  
$EndComp
Text GLabel 3125 4875 0    10   UnSpc ~ 0
GND
Text GLabel 3025 4875 2    10   UnSpc ~ 0
GND
Text GLabel 3125 4975 0    10   UnSpc ~ 0
GND
Text GLabel 3025 4975 2    10   UnSpc ~ 0
GND
Text GLabel 5625 5575 2    10   UnSpc ~ 0
GND
Text GLabel 5725 5575 0    10   UnSpc ~ 0
GND
Text GLabel 5625 5375 2    10   UnSpc ~ 0
GND
Text GLabel 5725 5375 0    10   UnSpc ~ 0
GND
Text GLabel 3725 5875 3    10   UnSpc ~ 0
GND
Text GLabel 3725 5975 1    10   UnSpc ~ 0
GND
Text GLabel 1200 6400 3    10   UnSpc ~ 0
GND
Text GLabel 1200 6500 1    10   UnSpc ~ 0
GND
Text GLabel 5125 6275 3    10   UnSpc ~ 0
GND
Text GLabel 5125 6375 1    10   UnSpc ~ 0
GND
Text GLabel 5375 7075 1    10   UnSpc ~ 0
GND
Text GLabel 5375 6975 3    10   UnSpc ~ 0
GND
Text GLabel 5625 4775 2    10   UnSpc ~ 0
GND
Text GLabel 5725 4775 0    10   UnSpc ~ 0
GND
Text Label 3125 4675 2    10   ~ 0
N$1
Text Label 3125 3875 2    10   ~ 0
N$2
Text Label 3125 3775 2    10   ~ 0
N$3
Text Label 3125 3975 2    10   ~ 0
N$4
Text Label 3125 4075 2    10   ~ 0
N$5
Text Label 3125 4175 2    10   ~ 0
N$6
Text Label 3125 4275 2    10   ~ 0
N$7
Text Label 3125 4375 2    10   ~ 0
N$8
Text Label 3125 4475 2    10   ~ 0
N$9
Text Label 3925 5875 3    10   ~ 0
N$11
Text Label 4325 3375 1    10   ~ 0
N$19
Text Label 4225 3375 1    10   ~ 0
N$20
Text Label 4125 3375 1    10   ~ 0
N$21
Text Label 3825 3375 1    10   ~ 0
N$23
Text Label 3925 3375 1    10   ~ 0
N$24
Text Label 4025 3375 1    10   ~ 0
N$25
Text Label 4425 3375 1    10   ~ 0
N$26
Text Label 4525 3375 1    10   ~ 0
N$27
Text Label 4625 3375 1    10   ~ 0
N$28
Text Label 4725 3375 1    10   ~ 0
N$29
Text Label 4825 3375 1    10   ~ 0
N$30
Text Label 4925 3375 1    10   ~ 0
N$31
Text Label 5025 3375 1    10   ~ 0
N$32
Text Label 3125 3675 2    10   ~ 0
N$33
Text Label 4025 5875 3    10   ~ 0
N$35
Text Label 4125 5875 3    10   ~ 0
N$36
Text Label 4225 5875 3    10   ~ 0
N$37
Text Label 4325 5875 3    10   ~ 0
N$38
Text Label 4425 5875 3    10   ~ 0
N$39
Text Label 4525 5875 3    10   ~ 0
N$40
Text Label 4625 5875 3    10   ~ 0
N$41
Text Label 5125 5875 3    10   ~ 0
N$42
Text Label 5125 5975 1    10   ~ 0
N$42
Text Label 5025 5875 3    10   ~ 0
N$43
Text Label 3125 5275 2    10   ~ 0
N$44
Text Label 3425 5875 3    10   ~ 0
N$46
Text Label 3525 5875 3    10   ~ 0
N$48
Text GLabel 5625 4975 2    10   UnSpc ~ 0
VCC
Text GLabel 5725 4975 0    10   UnSpc ~ 0
VCC
Text GLabel 4825 5875 3    10   UnSpc ~ 0
VCC
Text GLabel 4825 5975 1    10   UnSpc ~ 0
VCC
Text GLabel 3125 5075 0    10   UnSpc ~ 0
VCC
Text GLabel 3025 5075 2    10   UnSpc ~ 0
VCC
Text GLabel 1200 6100 1    10   UnSpc ~ 0
VCC
Text GLabel 1200 6000 3    10   UnSpc ~ 0
VCC
Text GLabel 5275 7075 1    10   UnSpc ~ 0
VCC
Text GLabel 5275 6975 3    10   UnSpc ~ 0
VCC
$Comp
L RHD2132_headstage-rescue:C_0402 C3
U 1 1 576BDBCF
P 1750 6225
F 0 "C3" H 1810 6239 70  0000 L BNN
F 1 "100n" H 1809 6040 70  0000 L BNN
F 2 "Resistors_SMD:R_0603" H 1740 6015 65  0001 L TNN
F 3 "" H 1750 6225 60  0001 C CNN
	1    1750 6225
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+015
U 1 1 576BDC25
P 1750 5900
F 0 "#P+015" V 1650 5800 70  0001 L BNN
F 1 "VCC" V 1650 5800 70  0000 L BNN
F 2 "" H 1750 5900 60  0001 C CNN
F 3 "" H 1750 5900 60  0001 C CNN
	1    1750 5900
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND016
U 1 1 576BDD07
P 1750 6675
F 0 "#GND016" H 1650 6575 70  0001 L BNN
F 1 "GND" H 1650 6575 70  0000 L BNN
F 2 "" H 1750 6675 60  0001 C CNN
F 3 "" H 1750 6675 60  0001 C CNN
	1    1750 6675
	1    0    0    -1  
$EndComp
Text Notes 3700 4275 0    118  Italic 24
RHD2132 Intan \n     Chip
Text Notes 3850 975  0    118  Italic 24
FFC Connectors
$Comp
L RHD2132_headstage-rescue:RHD2132_QFN U1
U 1 1 5AC64383
P 4325 4675
F 0 "U1" H 4125 4925 45  0000 L BNN
F 1 "RHD2132_QFN" H 4125 4725 45  0000 L BNN
F 2 "NNClib:RHD2000-QFN56_8_X_8_LARGE" H 4355 4825 20  0001 C CNN
F 3 "" H 4325 4675 60  0000 C CNN
	1    4325 4675
	1    0    0    -1  
$EndComp
Text GLabel 2925 4675 0    39   Input ~ 0
REF_ELECTRODE
$Comp
L RHD2132_headstage-rescue:FFCfemale20 J1
U 1 1 5AC65309
P 2900 1375
F 0 "J1" H 3150 1375 50  0000 C CNN
F 1 "FFCfemale20" H 3200 175 50  0000 C CNN
F 2 "NNClib:FFCfemale20" H 3100 775 50  0001 C CNN
F 3 "" H 3100 775 50  0001 C CNN
	1    2900 1375
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:FFCfemale20 J2
U 1 1 5AC653A2
P 5125 1375
F 0 "J2" H 5375 1375 50  0000 C CNN
F 1 "FFCfemale20" H 5425 175 50  0000 C CNN
F 2 "NNClib:FFCfemale20" H 5325 775 50  0001 C CNN
F 3 "" H 5325 775 50  0001 C CNN
	1    5125 1375
	1    0    0    -1  
$EndComp
$Comp
L RHD2132_headstage-rescue:Conn_01x08 J4
U 1 1 5AC65C1A
P 4250 6575
F 0 "J4" H 4250 6975 50  0000 C CNN
F 1 "SPI" H 4250 6075 50  0000 C CNN
F 2 "NNClib:JumperSPI2" H 4250 6575 50  0001 C CNN
F 3 "" H 4250 6575 50  0001 C CNN
	1    4250 6575
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:Conn_01x02 J5
U 1 1 5AC65C61
P 5375 7325
F 0 "J5" H 5375 7425 50  0000 C CNN
F 1 "POWER" H 5375 7125 50  0000 C CNN
F 2 "Pin_Headers:Pin_Header_Straight_1x02_Pitch2.00mm" H 5375 7325 50  0001 C CNN
F 3 "" H 5375 7325 50  0001 C CNN
	1    5375 7325
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:Conn_01x03 J3
U 1 1 5AC65EAC
P 3325 7300
F 0 "J3" H 3325 7500 50  0000 C CNN
F 1 "AuxInputs" H 3325 7100 50  0000 C CNN
F 2 "Pin_Headers:Pin_Header_Straight_1x03_Pitch2.00mm" H 3325 7300 50  0001 C CNN
F 3 "" H 3325 7300 50  0001 C CNN
	1    3325 7300
	0    1    1    0   
$EndComp
$Comp
L Device:R R1
U 1 1 5AC66076
P 6425 5100
F 0 "R1" V 6505 5100 50  0000 C CNN
F 1 "0" V 6425 5100 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 6355 5100 50  0001 C CNN
F 3 "" H 6425 5100 50  0001 C CNN
	1    6425 5100
	0    1    1    0   
$EndComp
$Comp
L Device:R R2
U 1 1 5AC660C3
P 6425 5300
F 0 "R2" V 6505 5300 50  0000 C CNN
F 1 "0" V 6425 5300 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 6355 5300 50  0001 C CNN
F 3 "" H 6425 5300 50  0001 C CNN
	1    6425 5300
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND017
U 1 1 5AC662CF
P 6775 5100
F 0 "#GND017" H 6675 5000 70  0001 L BNN
F 1 "GND" V 6800 4875 70  0000 L BNN
F 2 "" H 6775 5100 60  0001 C CNN
F 3 "" H 6775 5100 60  0001 C CNN
	1    6775 5100
	0    -1   -1   0   
$EndComp
$Comp
L RHD2132_headstage-rescue:VCC-RESCUE-RHD2132_headstage1 #P+018
U 1 1 5AC663A7
P 6775 5300
F 0 "#P+018" V 6675 5200 70  0001 L BNN
F 1 "VCC" V 6800 5325 70  0000 L BNN
F 2 "" H 6775 5300 60  0001 C CNN
F 3 "" H 6775 5300 60  0001 C CNN
	1    6775 5300
	0    1    -1   0   
$EndComp
Text Notes 6200 5000 0    39   Italic 8
Select between\nLVDS and CMOS
Text GLabel 3050 4475 0    39   Input ~ 0
CH0
Text GLabel 3050 4375 0    39   Input ~ 0
CH1
Text GLabel 3050 4275 0    39   Input ~ 0
CH2
Text GLabel 3050 4175 0    39   Input ~ 0
CH3
Text GLabel 3050 4075 0    39   Input ~ 0
CH4
Text GLabel 3050 3975 0    39   Input ~ 0
CH5
Text GLabel 3050 3875 0    39   Input ~ 0
CH6
Text GLabel 3050 3775 0    39   Input ~ 0
CH7
Text GLabel 3050 3675 0    39   Input ~ 0
CH8
Text GLabel 3725 3275 1    39   Input ~ 0
CH9
Text GLabel 3825 3300 1    39   Input ~ 0
CH10
Text GLabel 3925 3300 1    39   Input ~ 0
CH11
Text GLabel 4025 3300 1    39   Input ~ 0
CH12
Text GLabel 4125 3300 1    39   Input ~ 0
CH13
Text GLabel 4225 3300 1    39   Input ~ 0
CH14
Text GLabel 4325 3300 1    39   Input ~ 0
CH15
Text GLabel 4425 3300 1    39   Input ~ 0
CH16
Text GLabel 4525 3300 1    39   Input ~ 0
CH17
Text GLabel 4625 3300 1    39   Input ~ 0
CH18
Text GLabel 4725 3300 1    39   Input ~ 0
CH19
Text GLabel 4825 3300 1    39   Input ~ 0
CH20
Text GLabel 4925 3300 1    39   Input ~ 0
CH21
Text GLabel 5025 3300 1    39   Input ~ 0
CH22
Text Label 3725 3375 1    10   ~ 0
N$22
Text GLabel 5725 3675 2    39   Input ~ 0
CH23
Text GLabel 5725 3775 2    39   Input ~ 0
CH24
Text GLabel 5725 3875 2    39   Input ~ 0
CH25
Text GLabel 5725 3975 2    39   Input ~ 0
CH26
Text GLabel 5725 4075 2    39   Input ~ 0
CH27
Text GLabel 5725 4175 2    39   Input ~ 0
CH28
Text GLabel 5725 4275 2    39   Input ~ 0
CH29
Text GLabel 5725 4375 2    39   Input ~ 0
CH30
Text GLabel 5725 4475 2    39   Input ~ 0
CH31
Text GLabel 3425 5950 3    39   Input ~ 0
AUXIN2
Text GLabel 3525 5950 3    39   Input ~ 0
AUXIN3
Text GLabel 3025 5275 0    39   Input ~ 0
AUXIN1
Text GLabel 3925 5950 3    39   Input ~ 0
~CS~-
Text GLabel 4025 5950 3    39   Input ~ 0
~CS~+
Text GLabel 4125 5950 3    39   Input ~ 0
SCLK-
Text GLabel 4225 5950 3    39   Input ~ 0
SCLK+
Text GLabel 4325 5950 3    39   Input ~ 0
MOSI-
Text GLabel 4425 5950 3    39   Input ~ 0
MOSI+
Text GLabel 4525 5950 3    39   Input ~ 0
MISO-
Text GLabel 4625 5950 3    39   Input ~ 0
MISO+
Text Notes 3900 7475 0    118  Italic 24
Jumpers
Text GLabel 3950 6850 3    39   Input ~ 0
~CS~-
Text GLabel 4050 6850 3    39   Input ~ 0
~CS~+
Text GLabel 4150 6850 3    39   Input ~ 0
SCLK-
Text GLabel 4250 6850 3    39   Input ~ 0
SCLK+
Text GLabel 4350 6850 3    39   Input ~ 0
MOSI-
Text GLabel 4450 6850 3    39   Input ~ 0
MOSI+
Text GLabel 4550 6850 3    39   Input ~ 0
MISO-
Text GLabel 4650 6850 3    39   Input ~ 0
MISO+
Text GLabel 3325 7000 1    39   Input ~ 0
AUXIN2
Text GLabel 3225 7000 1    39   Input ~ 0
AUXIN3
Text GLabel 3425 7000 1    39   Input ~ 0
AUXIN1
Text GLabel 2825 2275 0    39   Input ~ 0
CH11
Text GLabel 2825 2175 0    39   Input ~ 0
CH10
Text GLabel 2825 2075 0    39   Input ~ 0
CH9
Text GLabel 2825 1975 0    39   Input ~ 0
CH8
Text GLabel 2825 1875 0    39   Input ~ 0
CH6
Text GLabel 2825 1775 0    39   Input ~ 0
CH4
Text GLabel 2825 1675 0    39   Input ~ 0
CH2
Text GLabel 2825 1575 0    39   Input ~ 0
CH0
Text GLabel 3475 1675 2    39   Input ~ 0
CH3
Text GLabel 3475 1775 2    39   Input ~ 0
CH5
Text GLabel 3475 1875 2    39   Input ~ 0
CH7
Text GLabel 3475 1975 2    39   Input ~ 0
CH15
Text GLabel 3475 2075 2    39   Input ~ 0
CH14
Text GLabel 3475 2175 2    39   Input ~ 0
CH12
Text GLabel 3475 2275 2    39   Input ~ 0
CH13
Text GLabel 5725 1575 2    39   Input ~ 0
CH30
Text GLabel 5725 1675 2    39   Input ~ 0
CH28
Text GLabel 5725 1775 2    39   Input ~ 0
CH26
Text GLabel 5725 1875 2    39   Input ~ 0
CH24
Text GLabel 5725 1975 2    39   Input ~ 0
CH22
Text GLabel 5725 2075 2    39   Input ~ 0
CH21
Text GLabel 5725 2175 2    39   Input ~ 0
CH20
Text GLabel 5725 2275 2    39   Input ~ 0
CH19
Text GLabel 5025 2275 0    39   Input ~ 0
CH17
Text GLabel 5025 2175 0    39   Input ~ 0
CH18
Text GLabel 5025 2075 0    39   Input ~ 0
CH16
Text GLabel 5025 1975 0    39   Input ~ 0
CH23
Text GLabel 5025 1875 0    39   Input ~ 0
CH25
Text GLabel 5025 1775 0    39   Input ~ 0
CH27
Text GLabel 5025 1675 0    39   Input ~ 0
CH29
Text GLabel 5025 1575 0    39   Input ~ 0
CH31
Wire Wire Line
	3950 6775 3950 6850
Wire Wire Line
	4050 6775 4050 6850
Wire Wire Line
	4150 6775 4150 6850
Wire Wire Line
	4250 6775 4250 6850
Wire Wire Line
	4350 6775 4350 6850
Wire Wire Line
	4450 6775 4450 6850
Wire Wire Line
	4550 6775 4550 6850
Wire Wire Line
	4650 6775 4650 6850
Wire Wire Line
	3425 7100 3425 7000
Wire Wire Line
	3325 7100 3325 7000
Wire Wire Line
	3225 7100 3225 7000
Wire Wire Line
	5875 5300 6275 5300
Wire Wire Line
	5875 5100 6275 5100
Wire Wire Line
	4625 5950 4625 5875
Wire Wire Line
	4525 5950 4525 5875
Wire Wire Line
	4425 5950 4425 5875
Wire Wire Line
	4325 5950 4325 5875
Wire Wire Line
	4225 5950 4225 5875
Wire Wire Line
	4125 5950 4125 5875
Wire Wire Line
	4025 5950 4025 5875
Wire Wire Line
	3925 5950 3925 5875
Wire Wire Line
	3025 5275 3125 5275
Wire Wire Line
	3525 5950 3525 5875
Wire Wire Line
	3425 5950 3425 5875
Wire Wire Line
	5725 4475 5625 4475
Wire Wire Line
	5725 4375 5625 4375
Wire Wire Line
	5725 4275 5625 4275
Wire Wire Line
	5725 4175 5625 4175
Wire Wire Line
	5725 4075 5625 4075
Wire Wire Line
	5725 3975 5625 3975
Wire Wire Line
	5725 3875 5625 3875
Wire Wire Line
	5725 3775 5625 3775
Wire Wire Line
	5725 3675 5625 3675
Wire Wire Line
	5025 3300 5025 3375
Wire Wire Line
	4925 3300 4925 3375
Wire Wire Line
	4825 3300 4825 3375
Wire Wire Line
	4725 3300 4725 3375
Wire Wire Line
	4625 3300 4625 3375
Wire Wire Line
	4525 3300 4525 3375
Wire Wire Line
	4425 3300 4425 3375
Wire Wire Line
	4325 3300 4325 3375
Wire Wire Line
	4225 3300 4225 3375
Wire Wire Line
	4125 3300 4125 3375
Wire Wire Line
	4025 3300 4025 3375
Wire Wire Line
	3925 3300 3925 3375
Wire Wire Line
	3825 3300 3825 3375
Wire Wire Line
	3725 3275 3725 3375
Wire Wire Line
	3050 3675 3125 3675
Wire Wire Line
	3050 3775 3125 3775
Wire Wire Line
	3050 3875 3125 3875
Wire Wire Line
	3050 3975 3125 3975
Wire Wire Line
	3050 4075 3125 4075
Wire Wire Line
	3050 4175 3125 4175
Wire Wire Line
	3050 4275 3125 4275
Wire Wire Line
	3050 4375 3125 4375
Wire Wire Line
	3050 4475 3125 4475
Wire Wire Line
	6575 5100 6675 5100
Wire Wire Line
	6575 5300 6675 5300
Connection ~ 5875 5175
Wire Wire Line
	5875 5100 5875 5175
Wire Wire Line
	2925 4675 3125 4675
Wire Wire Line
	1750 6425 1750 6575
Wire Wire Line
	1750 6000 1750 6125
Wire Wire Line
	5275 6975 5275 7125
Wire Wire Line
	1200 6000 1200 6100
Wire Wire Line
	3025 5075 3125 5075
Wire Wire Line
	4825 5975 4825 5875
Wire Wire Line
	5625 5175 5875 5175
Wire Wire Line
	5725 4975 5625 4975
Wire Wire Line
	5125 5975 5125 5875
Wire Wire Line
	5725 4775 5625 4775
Wire Wire Line
	5375 6975 5375 7125
Wire Wire Line
	5125 6375 5125 6275
Wire Wire Line
	1200 6500 1200 6400
Wire Wire Line
	3725 5975 3725 5875
Wire Wire Line
	5725 5375 5625 5375
Wire Wire Line
	5725 5575 5625 5575
Wire Wire Line
	3025 4975 3125 4975
Wire Wire Line
	3025 4875 3125 4875
Wire Wire Line
	2825 1575 2900 1575
Wire Wire Line
	2825 1675 2900 1675
Wire Wire Line
	2825 1775 2900 1775
Wire Wire Line
	2825 1875 2900 1875
Wire Wire Line
	2825 1975 2900 1975
Wire Wire Line
	2825 2075 2900 2075
Wire Wire Line
	2825 2175 2900 2175
Wire Wire Line
	2825 2275 2900 2275
Wire Wire Line
	3475 2275 3400 2275
Wire Wire Line
	3475 2175 3400 2175
Wire Wire Line
	3475 2075 3400 2075
Wire Wire Line
	3475 1975 3400 1975
Wire Wire Line
	3400 1875 3475 1875
Wire Wire Line
	3400 1775 3475 1775
Wire Wire Line
	3400 1675 3475 1675
Wire Wire Line
	5025 1575 5125 1575
Wire Wire Line
	5025 1675 5125 1675
Wire Wire Line
	5025 1775 5125 1775
Wire Wire Line
	5025 1875 5125 1875
Wire Wire Line
	5025 1975 5125 1975
Wire Wire Line
	5025 2075 5125 2075
Wire Wire Line
	5025 2175 5125 2175
Wire Wire Line
	5025 2275 5125 2275
Wire Wire Line
	5725 1575 5625 1575
Wire Wire Line
	5725 1675 5625 1675
Wire Wire Line
	5725 1775 5625 1775
Wire Wire Line
	5725 1875 5625 1875
Wire Wire Line
	5725 1975 5625 1975
Wire Wire Line
	5725 2075 5625 2075
Wire Wire Line
	5725 2175 5625 2175
Wire Wire Line
	5725 2275 5625 2275
Wire Wire Line
	3400 1575 3475 1575
Text GLabel 3475 1575 2    39   Input ~ 0
CH1
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND019
U 1 1 5AD3B9E7
P 2700 1425
F 0 "#GND019" H 2601 1325 70  0001 L BNN
F 1 "GND" V 2725 1200 70  0000 L BNN
F 2 "" H 2700 1425 60  0001 C CNN
F 3 "" H 2700 1425 60  0001 C CNN
	1    2700 1425
	0    1    1    0   
$EndComp
Text GLabel 3475 1475 2    39   Input ~ 0
REF_ELECTRODE
Wire Wire Line
	2800 1475 2900 1475
Wire Wire Line
	3475 1475 3400 1475
Text GLabel 2825 2375 0    39   Input ~ 0
REF_ELECTRODE
Wire Wire Line
	3400 2375 3525 2375
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND020
U 1 1 5AD3C0DC
P 3625 2425
F 0 "#GND020" H 3526 2325 70  0001 L BNN
F 1 "GND" V 3650 2200 70  0000 L BNN
F 2 "" H 3625 2425 60  0001 C CNN
F 3 "" H 3625 2425 60  0001 C CNN
	1    3625 2425
	0    -1   -1   0   
$EndComp
Text GLabel 5725 1475 2    39   Input ~ 0
REF_ELECTRODE
Wire Wire Line
	5025 1475 5125 1475
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND021
U 1 1 5AD3C355
P 4925 1425
F 0 "#GND021" H 4826 1325 70  0001 L BNN
F 1 "GND" V 4950 1200 70  0000 L BNN
F 2 "" H 4925 1425 60  0001 C CNN
F 3 "" H 4925 1425 60  0001 C CNN
	1    4925 1425
	0    1    1    0   
$EndComp
Wire Wire Line
	5725 1475 5625 1475
Text GLabel 5025 2375 0    39   Input ~ 0
REF_ELECTRODE
Wire Wire Line
	5725 2375 5625 2375
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND022
U 1 1 5AD3C546
P 5825 2425
F 0 "#GND022" H 5726 2325 70  0001 L BNN
F 1 "GND" V 5850 2200 70  0000 L BNN
F 2 "" H 5825 2425 60  0001 C CNN
F 3 "" H 5825 2425 60  0001 C CNN
	1    5825 2425
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5025 2375 5125 2375
Wire Wire Line
	2800 1425 2800 1475
Wire Wire Line
	2900 2375 2825 2375
Wire Wire Line
	3525 2375 3525 2425
Wire Wire Line
	5025 1425 5025 1475
Wire Wire Line
	5725 2425 5725 2375
$Comp
L RHD2132_headstage-rescue:Conn_01x01 J6
U 1 1 5AE328E4
P 7725 4325
F 0 "J6" H 7725 4425 50  0000 C CNN
F 1 "Conn_01x01" H 7725 4225 50  0000 C CNN
F 2 "Pin_Headers:Pin_Header_Straight_1x01_Pitch2.00mm" H 7725 4325 50  0001 C CNN
F 3 "" H 7725 4325 50  0001 C CNN
	1    7725 4325
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:Conn_01x01 J7
U 1 1 5AE3299B
P 8300 4325
F 0 "J7" H 8300 4425 50  0000 C CNN
F 1 "Conn_01x01" H 8300 4225 50  0000 C CNN
F 2 "Pin_Headers:Pin_Header_Straight_1x01_Pitch2.00mm" H 8300 4325 50  0001 C CNN
F 3 "" H 8300 4325 50  0001 C CNN
	1    8300 4325
	0    1    1    0   
$EndComp
Text Notes 7225 3900 0    118  Italic 24
Mechanical support
Text Notes 7000 5375 0    39   Italic 8
LVDS
Text Notes 7000 5100 0    39   Italic 8
CMOS
Text GLabel 1600 4650 0    39   Input ~ 0
REF_ELECTRODE
$Comp
L Device:R R3
U 1 1 5B06A65F
P 1750 4650
F 0 "R3" V 1830 4650 50  0000 C CNN
F 1 "0" V 1750 4650 50  0000 C CNN
F 2 "Resistors_SMD:R_0603" V 1680 4650 50  0001 C CNN
F 3 "" H 1750 4650 50  0001 C CNN
	1    1750 4650
	0    1    1    0   
$EndComp
$Comp
L RHD2132_headstage-rescue:GND-RESCUE-RHD2132_headstage1 #GND023
U 1 1 5B06A79A
P 2000 4650
F 0 "#GND023" H 1901 4550 70  0001 L BNN
F 1 "GND" V 2000 4425 70  0000 L BNN
F 2 "" H 2000 4650 60  0001 C CNN
F 3 "" H 2000 4650 60  0001 C CNN
	1    2000 4650
	0    -1   -1   0   
$EndComp
Text Notes 1275 4500 0    39   Italic 8
    Select between\nRHD2216 and RHD2132
Text Notes 1250 4975 0    39   Italic 8
This resistor must be soldered\n    when using RHD2216.
Text Notes 5650 4675 0    39   Italic 8
Unconnected
Text Notes 5000 5900 3    39   Italic 8
Unconnected
Wire Wire Line
	5875 5175 5875 5300
$EndSCHEMATC
