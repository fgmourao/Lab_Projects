EESchema Schematic File Version 4
LIBS:NNC_monitor-cache
EELAYER 26 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L Amplifier_Operational:LM324 U1
U 1 1 5C4DC812
P 3650 3200
F 0 "U1" H 3650 2833 50  0000 C CNN
F 1 "LM324" H 3650 2924 50  0000 C CNN
F 2 "Housings_DIP:DIP-14_W7.62mm" H 3600 3300 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/lm2902-n.pdf" H 3700 3400 50  0001 C CNN
	1    3650 3200
	1    0    0    1   
$EndComp
$Comp
L Device:R R3
U 1 1 5C4DF616
P 2500 3000
F 0 "R3" H 2570 3046 50  0000 L CNN
F 1 "100K" H 2570 2955 50  0000 L CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 2430 3000 50  0001 C CNN
F 3 "~" H 2500 3000 50  0001 C CNN
	1    2500 3000
	1    0    0    -1  
$EndComp
$Comp
L Device:R R4
U 1 1 5C4DFBBC
P 2500 3650
F 0 "R4" H 2570 3696 50  0000 L CNN
F 1 "100K" H 2570 3605 50  0000 L CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 2430 3650 50  0001 C CNN
F 3 "~" H 2500 3650 50  0001 C CNN
	1    2500 3650
	1    0    0    -1  
$EndComp
Wire Wire Line
	3350 3300 2500 3300
Wire Wire Line
	2500 3150 2500 3300
Wire Wire Line
	2500 3300 2500 3500
Connection ~ 2500 3300
Wire Wire Line
	3950 3200 4050 3200
Wire Wire Line
	4050 3200 4050 2650
Wire Wire Line
	4050 2650 3250 2650
Wire Wire Line
	3250 2650 3250 3100
Wire Wire Line
	3250 3100 3350 3100
$Comp
L Amplifier_Operational:LM324 U1
U 2 1 5C4E2B40
P 3650 2150
F 0 "U1" H 3650 1783 50  0000 C CNN
F 1 "LM324" H 3650 1874 50  0000 C CNN
F 2 "Housings_DIP:DIP-14_W7.62mm" H 3600 2250 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/lm2902-n.pdf" H 3700 2350 50  0001 C CNN
	2    3650 2150
	1    0    0    1   
$EndComp
$Comp
L Device:R R5
U 1 1 5C4E4C73
P 3650 1600
F 0 "R5" V 3443 1600 50  0000 C CNN
F 1 "82k" V 3534 1600 50  0000 C CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 3580 1600 50  0001 C CNN
F 3 "~" H 3650 1600 50  0001 C CNN
	1    3650 1600
	0    1    1    0   
$EndComp
$Comp
L Device:C C1
U 1 1 5C4E5582
P 3650 1200
F 0 "C1" V 3398 1200 50  0000 C CNN
F 1 "100nF" V 3489 1200 50  0000 C CNN
F 2 "Capacitors_ThroughHole:C_Disc_D3_P2.5" H 3688 1050 50  0001 C CNN
F 3 "~" H 3650 1200 50  0001 C CNN
	1    3650 1200
	0    1    1    0   
$EndComp
Wire Wire Line
	3350 2250 3250 2250
Wire Wire Line
	3250 2250 3250 2450
Connection ~ 3250 2650
Wire Wire Line
	3800 1200 4050 1200
Wire Wire Line
	4050 2150 3950 2150
Wire Wire Line
	3800 1600 4050 1600
Connection ~ 4050 1600
Wire Wire Line
	4050 1600 4050 2150
Wire Wire Line
	3250 2050 3250 1600
Wire Wire Line
	3250 1200 3500 1200
Wire Wire Line
	3250 1600 3500 1600
Connection ~ 3250 1600
Wire Wire Line
	3250 1600 3250 1200
$Comp
L Amplifier_Operational:LM324 U1
U 3 1 5C4E8E2B
P 5200 2150
F 0 "U1" H 5200 1783 50  0000 C CNN
F 1 "LM324" H 5200 1874 50  0000 C CNN
F 2 "Housings_DIP:DIP-14_W7.62mm" H 5150 2250 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/lm2902-n.pdf" H 5250 2350 50  0001 C CNN
	3    5200 2150
	1    0    0    1   
$EndComp
Wire Wire Line
	3250 2450 4800 2450
Wire Wire Line
	4800 2450 4800 2250
Wire Wire Line
	4800 2250 4900 2250
Connection ~ 3250 2450
Wire Wire Line
	3250 2450 3250 2650
$Comp
L Device:R_POT_TRIM RV1
U 1 1 5C4EF563
P 6350 1550
F 0 "RV1" H 6280 1596 50  0000 R CNN
F 1 "1k-10k" H 6280 1505 50  0000 R CNN
F 2 "Potentiometers:Potentiometer_Trimmer-Suntan-TSR-3386P" H 6350 1550 50  0001 C CNN
F 3 "~" H 6350 1550 50  0001 C CNN
	1    6350 1550
	-1   0    0    -1  
$EndComp
Wire Wire Line
	5500 2150 5600 2150
Wire Wire Line
	4900 2050 4900 1750
$Comp
L Amplifier_Operational:LM324 U1
U 5 1 5C4F38AE
P 1750 1000
F 0 "U1" H 1708 1046 50  0000 L CNN
F 1 "LM324" H 1708 955 50  0000 L CNN
F 2 "Housings_DIP:DIP-14_W7.62mm" H 1700 1100 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/lm2902-n.pdf" H 1800 1200 50  0001 C CNN
	5    1750 1000
	0    -1   -1   0   
$EndComp
$Comp
L Connector:Conn_01x04_Male J4
U 1 1 5C51766B
P 1800 6850
F 0 "J4" H 1772 6824 50  0000 R CNN
F 1 "OLED Monitor" H 1772 6733 50  0000 R CNN
F 2 "Connector_Molex:Molex_KK-254_AE-6410-04A_1x04_P2.54mm_Vertical" H 1800 6850 50  0001 C CNN
F 3 "~" H 1800 6850 50  0001 C CNN
	1    1800 6850
	-1   0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x03_Male J3
U 1 1 5C519CA3
P 1850 6000
F 0 "J3" H 1822 6024 50  0000 R CNN
F 1 "Relay" H 1822 5933 50  0000 R CNN
F 2 "Connector_Molex:Molex_KK-254_AE-6410-03A_1x03_P2.54mm_Vertical" H 1850 6000 50  0001 C CNN
F 3 "~" H 1850 6000 50  0001 C CNN
	1    1850 6000
	-1   0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J2
U 1 1 5C54C041
P 1850 5150
F 0 "J2" H 1822 5124 50  0000 R CNN
F 1 "Thermistor" H 1822 5033 50  0000 R CNN
F 2 "Connector_Molex:Molex_KK-254_AE-6410-02A_1x02_P2.54mm_Vertical" H 1850 5150 50  0001 C CNN
F 3 "~" H 1850 5150 50  0001 C CNN
	1    1850 5150
	-1   0    0    -1  
$EndComp
$Comp
L Device:R R6
U 1 1 5C550BD8
P 1650 4750
F 0 "R6" H 1720 4796 50  0000 L CNN
F 1 "10K" H 1720 4705 50  0000 L CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 1580 4750 50  0001 C CNN
F 3 "~" H 1650 4750 50  0001 C CNN
	1    1650 4750
	1    0    0    -1  
$EndComp
Wire Wire Line
	1650 4400 1650 4600
$Comp
L Device:C C2
U 1 1 5C4E7D31
P 5200 1550
F 0 "C2" V 4948 1550 50  0000 C CNN
F 1 "100nF" V 5039 1550 50  0000 C CNN
F 2 "Capacitors_ThroughHole:C_Disc_D3_P2.5" H 5238 1400 50  0001 C CNN
F 3 "~" H 5200 1550 50  0001 C CNN
	1    5200 1550
	0    1    1    0   
$EndComp
Wire Wire Line
	4900 1550 5050 1550
Wire Wire Line
	5350 1550 5600 1550
Wire Wire Line
	5600 1550 5600 2150
Connection ~ 5600 2150
Wire Wire Line
	5600 2150 6350 2150
Wire Wire Line
	4050 1200 4050 1300
Wire Wire Line
	6350 1400 6350 900 
Wire Wire Line
	6350 900  4250 900 
Wire Wire Line
	4250 900  4250 1300
Wire Wire Line
	4250 1300 4050 1300
Connection ~ 4050 1300
Wire Wire Line
	4050 1300 4050 1600
Wire Wire Line
	4900 1750 4700 1750
Wire Wire Line
	4700 1750 4700 1050
Wire Wire Line
	4700 1050 5900 1050
Wire Wire Line
	5900 1050 5900 1550
Wire Wire Line
	5900 1550 6200 1550
Connection ~ 4900 1750
Wire Wire Line
	4900 1750 4900 1550
$Comp
L Device:CP1 C3
U 1 1 5C4F946F
P 2650 2050
F 0 "C3" V 2902 2050 50  0000 C CNN
F 1 "CP1" V 2811 2050 50  0000 C CNN
F 2 "Capacitors_ThroughHole:C_Radial_D5_L6_P2.5" H 2650 2050 50  0001 C CNN
F 3 "~" H 2650 2050 50  0001 C CNN
	1    2650 2050
	0    -1   -1   0   
$EndComp
$Comp
L Connector:Conn_01x04_Male J1
U 1 1 5C500CBB
P 800 2150
F 0 "J1" H 908 2431 50  0000 C CNN
F 1 "Optical Sensor" H 908 2340 50  0000 C CNN
F 2 "Connector_Molex:Molex_KK-254_AE-6410-04A_1x04_P2.54mm_Vertical" H 800 2150 50  0001 C CNN
F 3 "~" H 800 2150 50  0001 C CNN
	1    800  2150
	1    0    0    -1  
$EndComp
$Comp
L Device:R R1
U 1 1 5C4F3E9A
P 1650 2250
F 0 "R1" H 1720 2296 50  0000 L CNN
F 1 "10k" H 1720 2205 50  0000 L CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 1580 2250 50  0001 C CNN
F 3 "~" H 1650 2250 50  0001 C CNN
	1    1650 2250
	0    1    1    0   
$EndComp
$Comp
L Device:R R2
U 1 1 5C4F457B
P 1650 2500
F 0 "R2" H 1720 2546 50  0000 L CNN
F 1 "220" H 1720 2455 50  0000 L CNN
F 2 "Resistors_ThroughHole:Resistor_Horizontal_RM10mm" V 1580 2500 50  0001 C CNN
F 3 "~" H 1650 2500 50  0001 C CNN
	1    1650 2500
	0    1    1    0   
$EndComp
Wire Wire Line
	3250 2050 3350 2050
Wire Wire Line
	2800 2050 3250 2050
Connection ~ 3250 2050
Wire Wire Line
	1000 2050 1350 2050
Wire Wire Line
	1800 2250 1900 2250
Wire Wire Line
	1800 2500 1900 2500
Wire Wire Line
	1500 2250 1350 2250
Wire Wire Line
	1350 2250 1350 2050
Connection ~ 1350 2050
Wire Wire Line
	1350 2050 2500 2050
Wire Wire Line
	1000 2250 1100 2250
Wire Wire Line
	1200 2150 1200 2500
Wire Wire Line
	1200 2500 1500 2500
Wire Wire Line
	1000 2150 1200 2150
NoConn ~ 15050 -3800
NoConn ~ -250 -950
Text GLabel 4400 4650 2    50   Output ~ 0
RESET
Wire Wire Line
	4250 4650 4400 4650
Text GLabel 4400 4850 2    50   Output ~ 0
IOREF
Text GLabel 4400 5050 2    50   Output ~ 0
AREF
Text GLabel 4400 5250 2    50   Input ~ 0
A0
Text GLabel 4400 5350 2    50   Input ~ 0
A1
Text GLabel 4400 5450 2    50   Input ~ 0
A2
Text GLabel 4400 5550 2    50   Input ~ 0
A3
Text GLabel 4400 5650 2    50   BiDi ~ 0
A4_SDA
Text GLabel 4400 5750 2    50   BiDi ~ 0
A5_SCL
Text GLabel 4400 5950 2    50   Input ~ 0
SDA
Text GLabel 4400 6050 2    50   Input ~ 0
SCL
Text GLabel 10450 18550 0    50   Input ~ 0
RESET
Text GLabel 10600 18550 0    50   Input ~ 0
RESET
Text GLabel -400 15200 0    50   Input ~ 0
GND
Text GLabel 3150 5950 0    50   Output ~ 0
D13
Text GLabel 3150 5850 0    50   Output ~ 0
D12
Text GLabel 3150 5750 0    50   Output ~ 0
D11
Text GLabel 3150 5650 0    50   Output ~ 0
D10
Text GLabel 3150 5550 0    50   Output ~ 0
D9
Text GLabel 3150 5450 0    50   Output ~ 0
D8
Text GLabel 3150 5350 0    50   Output ~ 0
D7
Text GLabel 3150 5250 0    50   Output ~ 0
D6
Text GLabel 3150 5150 0    50   Output ~ 0
D5
Text GLabel 3150 5050 0    50   Output ~ 0
D4
Text GLabel 3150 4950 0    50   Output ~ 0
D3
Text GLabel 3150 4850 0    50   Output ~ 0
D2
Text GLabel 3150 4750 0    50   BiDi ~ 0
D1_TX
Text GLabel 3150 4650 0    50   BiDi ~ 0
D0_RX
Wire Wire Line
	3150 4650 3250 4650
Wire Wire Line
	3150 4750 3250 4750
Wire Wire Line
	3150 5050 3250 5050
Wire Wire Line
	3150 5150 3250 5150
Wire Wire Line
	3150 5250 3250 5250
Wire Wire Line
	3150 5350 3250 5350
Wire Wire Line
	3150 5450 3250 5450
Wire Wire Line
	3150 5550 3250 5550
Wire Wire Line
	3150 5650 3250 5650
Wire Wire Line
	3150 5750 3250 5750
Wire Wire Line
	3150 5850 3250 5850
Wire Wire Line
	3150 5950 3250 5950
Wire Wire Line
	4250 4850 4400 4850
Wire Wire Line
	4250 5050 4400 5050
Wire Wire Line
	4250 5450 4400 5450
Wire Wire Line
	4250 5550 4400 5550
Wire Wire Line
	4250 5950 4400 5950
Wire Wire Line
	4250 6050 4400 6050
Text GLabel 3850 4150 1    50   Output ~ 0
3V3
Text GLabel 3950 4150 1    50   Output ~ 0
+5V
Wire Wire Line
	3150 4950 3250 4950
Wire Wire Line
	1100 2250 1100 2750
Wire Wire Line
	1000 2350 1000 2750
Text GLabel 1900 2250 2    50   Input ~ 0
+5V
Text GLabel 1900 2500 2    50   Input ~ 0
+5V
Text GLabel 1100 2750 3    50   Input ~ 0
GND
Text GLabel 1000 2750 3    50   Input ~ 0
GND
Text GLabel 1300 1100 0    50   Input ~ 0
+5V
Text GLabel 2100 3800 0    50   Input ~ 0
GND
Wire Wire Line
	2100 3800 2500 3800
Wire Wire Line
	2100 2850 2500 2850
Text GLabel 2100 2850 0    50   Input ~ 0
+5V
Wire Wire Line
	1300 1100 1450 1100
Text GLabel 6650 1900 2    50   Output ~ 0
A0
Wire Wire Line
	6350 1700 6350 1900
Wire Wire Line
	6350 1900 6650 1900
Connection ~ 6350 1900
Wire Wire Line
	6350 1900 6350 2150
Wire Wire Line
	4250 5250 4400 5250
Wire Wire Line
	4250 5350 4400 5350
Text GLabel 1500 5150 0    50   Output ~ 0
A1
Wire Wire Line
	1500 5150 1650 5150
Text GLabel 1650 4400 1    50   Input ~ 0
+5V
Text GLabel 1500 5250 0    50   Input ~ 0
GND
Wire Wire Line
	1500 5250 1650 5250
Text GLabel 1500 5900 0    50   Input ~ 0
+5V
Wire Wire Line
	1500 5900 1650 5900
Text GLabel 1500 6000 0    50   Input ~ 0
D2
Wire Wire Line
	3150 4850 3250 4850
Wire Wire Line
	1500 6000 1650 6000
Text GLabel 1500 6100 0    50   Input ~ 0
GND
Wire Wire Line
	1500 6100 1650 6100
Text GLabel 1450 6850 0    50   Output ~ 0
A4_SDA
Wire Wire Line
	1450 6850 1600 6850
Text GLabel 1450 6950 0    50   Output ~ 0
A5_SCL
Wire Wire Line
	1450 6950 1600 6950
Wire Wire Line
	4250 5650 4400 5650
Wire Wire Line
	4250 5750 4400 5750
Text GLabel 1450 7050 0    50   Input ~ 0
GND
Wire Wire Line
	1450 7050 1600 7050
Text GLabel 1450 6750 0    50   Input ~ 0
+5V
Wire Wire Line
	1450 6750 1600 6750
Wire Wire Line
	1650 4900 1650 5150
Connection ~ 1650 5150
Text GLabel 6000 5500 0    50   Input ~ 0
D13
Text GLabel 6000 5400 0    50   Input ~ 0
D12
Text GLabel 6000 5300 0    50   Input ~ 0
D11
Text GLabel 6000 5200 0    50   Input ~ 0
D10
Text GLabel 6000 5100 0    50   Input ~ 0
D9
Text GLabel 6000 5000 0    50   Input ~ 0
D8
Text GLabel 6000 4850 0    50   Input ~ 0
D7
Text GLabel 6000 4750 0    50   Input ~ 0
D6
Text GLabel 6000 4650 0    50   Input ~ 0
D5
Text GLabel 6000 4550 0    50   Input ~ 0
D4
Text GLabel 6000 4450 0    50   Input ~ 0
D3
Text GLabel 6000 4350 0    50   Input ~ 0
D2
Text GLabel 6000 4250 0    50   BiDi ~ 0
D1_TX
Text GLabel 6000 4150 0    50   BiDi ~ 0
D0_RX
Wire Wire Line
	6000 4150 6200 4150
Wire Wire Line
	6000 4250 6200 4250
Wire Wire Line
	6000 4350 6200 4350
Wire Wire Line
	6000 4450 6200 4450
Wire Wire Line
	6000 4550 6200 4550
Wire Wire Line
	6000 4650 6200 4650
Wire Wire Line
	6000 4750 6200 4750
Wire Wire Line
	6000 4850 6200 4850
Wire Wire Line
	6000 5000 6200 5000
Wire Wire Line
	6000 5100 6200 5100
Wire Wire Line
	6000 5200 6200 5200
Wire Wire Line
	6000 5300 6200 5300
Wire Wire Line
	6000 5400 6200 5400
Wire Wire Line
	6000 5500 6200 5500
Text GLabel 6200 6400 2    50   Input ~ 0
RESET
Text GLabel 6200 6500 2    50   Input ~ 0
IOREF
Text GLabel 6250 5950 2    50   Input ~ 0
AREF
Text GLabel 6250 7050 2    50   Output ~ 0
A0
Text GLabel 6250 7150 2    50   Output ~ 0
A1
Text GLabel 6250 7250 2    50   Output ~ 0
A2
Text GLabel 6250 7350 2    50   Output ~ 0
A3
Text GLabel 6250 7450 2    50   BiDi ~ 0
A4_SDA
Text GLabel 6250 7550 2    50   BiDi ~ 0
A5_SCL
Text GLabel 6250 5850 2    50   Input ~ 0
SDA
Text GLabel 6250 5750 2    50   Input ~ 0
SCL
Wire Wire Line
	6000 5750 6250 5750
Wire Wire Line
	6000 5850 6250 5850
Wire Wire Line
	6000 5950 6250 5950
Wire Wire Line
	6000 7050 6250 7050
Wire Wire Line
	6000 7150 6250 7150
Wire Wire Line
	6000 7250 6250 7250
Wire Wire Line
	6000 7350 6250 7350
Wire Wire Line
	6000 7450 6250 7450
Wire Wire Line
	6000 7550 6250 7550
Wire Wire Line
	3850 4250 3850 4150
Wire Wire Line
	3950 4250 3950 4150
$Comp
L Connector:Conn_01x03_Female J7
U 1 1 5C9651B2
P 6350 3750
F 0 "J7" H 6378 3776 50  0000 L CNN
F 1 "Conn_01x03_Female" H 6378 3685 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x03_P2.54mm_Vertical" H 6350 3750 50  0001 C CNN
F 3 "~" H 6350 3750 50  0001 C CNN
	1    6350 3750
	1    0    0    -1  
$EndComp
Text GLabel 6000 3100 0    50   Input ~ 0
+5V
Text GLabel 6000 2750 0    50   Input ~ 0
3V3
Text GLabel 6000 3650 0    50   Input ~ 0
GND
Text GLabel 6000 3750 0    50   Input ~ 0
GND
Text GLabel 6000 3850 0    50   Input ~ 0
GND
Wire Wire Line
	6000 3650 6150 3650
Wire Wire Line
	6000 3750 6150 3750
Wire Wire Line
	6000 3850 6150 3850
NoConn ~ 3650 4250
$Comp
L Connector:Conn_01x02_Female J10
U 1 1 5C9B632C
P 6350 3100
F 0 "J10" H 6378 3076 50  0000 L CNN
F 1 "Conn_01x02_Female" H 6378 2985 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x02_P2.54mm_Vertical" H 6350 3100 50  0001 C CNN
F 3 "~" H 6350 3100 50  0001 C CNN
	1    6350 3100
	1    0    0    -1  
$EndComp
Wire Wire Line
	6000 3100 6150 3100
Wire Wire Line
	6150 2750 6000 2750
Text GLabel 6000 3200 0    50   Input ~ 0
+5V
Text GLabel 6000 2850 0    50   Input ~ 0
3V3
$Comp
L Connector:Conn_01x02_Female J9
U 1 1 5C9C5B4D
P 6350 2750
F 0 "J9" H 6378 2726 50  0000 L CNN
F 1 "Conn_01x02_Female" H 6378 2635 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x02_P2.54mm_Vertical" H 6350 2750 50  0001 C CNN
F 3 "~" H 6350 2750 50  0001 C CNN
	1    6350 2750
	1    0    0    -1  
$EndComp
Wire Wire Line
	6000 3200 6150 3200
Wire Wire Line
	6150 2850 6000 2850
Text GLabel 7450 3100 0    50   Input ~ 0
+5V
Text GLabel 7450 2700 0    50   Input ~ 0
3V3
$Comp
L Connector:Conn_01x02_Female J8
U 1 1 5C9CB3A6
P 7800 3100
F 0 "J8" H 7828 3076 50  0000 L CNN
F 1 "Conn_01x02_Female" H 7828 2985 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x02_P2.54mm_Vertical" H 7800 3100 50  0001 C CNN
F 3 "~" H 7800 3100 50  0001 C CNN
	1    7800 3100
	1    0    0    -1  
$EndComp
Wire Wire Line
	7450 3100 7600 3100
Wire Wire Line
	7600 2700 7450 2700
Wire Wire Line
	3650 6450 3650 6350
Text GLabel 3650 6450 3    50   Output ~ 0
GND
Wire Wire Line
	2050 1100 2200 1100
$Comp
L MCU_Module:Arduino_UNO_R3 A1
U 1 1 5C918044
P 3750 5250
F 0 "A1" H 3750 6431 50  0000 C CNN
F 1 "Arduino_UNO_R3" H 3750 6340 50  0000 C CNN
F 2 "Module:Arduino_UNO_R3" H 3900 4200 50  0001 L CNN
F 3 "https://www.arduino.cc/en/Main/arduinoBoardUno" H 3550 6300 50  0001 C CNN
	1    3750 5250
	1    0    0    -1  
$EndComp
NoConn ~ 3850 6350
Text GLabel 2200 1100 2    50   Input ~ 0
GND
NoConn ~ 3750 6350
$Comp
L Connector:Conn_01x03_Female J5
U 1 1 5C4FA30C
P 5800 5850
F 0 "J5" H 5828 5876 50  0000 L CNN
F 1 "Conn_01x03_Female" H 5828 5785 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x03_P2.54mm_Vertical" H 5800 5850 50  0001 C CNN
F 3 "~" H 5800 5850 50  0001 C CNN
	1    5800 5850
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J12
U 1 1 5C5421FF
P 5700 6500
F 0 "J12" H 5728 6476 50  0000 L CNN
F 1 "Conn_01x02_Female" H 5728 6385 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x02_P2.54mm_Vertical" H 5700 6500 50  0001 C CNN
F 3 "~" H 5700 6500 50  0001 C CNN
	1    5700 6500
	-1   0    0    1   
$EndComp
Wire Wire Line
	5900 6400 6200 6400
Wire Wire Line
	5900 6500 6200 6500
Text GLabel 7450 3200 0    50   Input ~ 0
+5V
Text GLabel 7450 2800 0    50   Input ~ 0
3V3
$Comp
L Connector:Conn_01x02_Female J13
U 1 1 5C53227C
P 7800 2700
F 0 "J13" H 7828 2676 50  0000 L CNN
F 1 "Conn_01x02_Female" H 7828 2585 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x02_P2.54mm_Vertical" H 7800 2700 50  0001 C CNN
F 3 "~" H 7800 2700 50  0001 C CNN
	1    7800 2700
	1    0    0    -1  
$EndComp
Wire Wire Line
	7450 3200 7600 3200
Wire Wire Line
	7600 2800 7450 2800
$Comp
L Connector:Conn_01x06_Female J11
U 1 1 5C54A5F0
P 5800 7350
F 0 "J11" H 5692 6825 50  0000 C CNN
F 1 "Conn_01x06_Female" H 5692 6916 50  0000 C CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x06_P2.54mm_Vertical" H 5800 7350 50  0001 C CNN
F 3 "~" H 5800 7350 50  0001 C CNN
	1    5800 7350
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x08_Female J6
U 1 1 5C587177
P 6400 4450
F 0 "J6" H 6428 4426 50  0000 L CNN
F 1 "Conn_01x08_Female" H 6428 4335 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x08_P2.54mm_Vertical" H 6400 4450 50  0001 C CNN
F 3 "~" H 6400 4450 50  0001 C CNN
	1    6400 4450
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x06_Female J14
U 1 1 5C589688
P 6400 5200
F 0 "J14" H 6428 5176 50  0000 L CNN
F 1 "Conn_01x06_Female" H 6428 5085 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x06_P2.54mm_Vertical" H 6400 5200 50  0001 C CNN
F 3 "~" H 6400 5200 50  0001 C CNN
	1    6400 5200
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x03_Female J15
U 1 1 5C5CB17E
P 8650 1500
F 0 "J15" H 8678 1526 50  0000 L CNN
F 1 "Conn_01x03_Female" H 8678 1435 50  0000 L CNN
F 2 "Connector_PinSocket_2.54mm:PinSocket_1x03_P2.54mm_Vertical" H 8650 1500 50  0001 C CNN
F 3 "~" H 8650 1500 50  0001 C CNN
	1    8650 1500
	1    0    0    -1  
$EndComp
Wire Wire Line
	7650 1600 7650 1850
Wire Wire Line
	7650 1850 8450 1850
Wire Wire Line
	8450 1850 8450 1600
Wire Wire Line
	7650 1400 7650 1050
Wire Wire Line
	7650 1050 8450 1050
Wire Wire Line
	8450 1050 8450 1400
Wire Wire Line
	8250 1500 8450 1500
$Comp
L Amplifier_Operational:LM324 U1
U 4 1 5C5F04FE
P 7950 1500
F 0 "U1" H 7950 1867 50  0000 C CNN
F 1 "LM324" H 7950 1776 50  0000 C CNN
F 2 "Housings_DIP:DIP-14_W7.62mm" H 7900 1600 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/lm2902-n.pdf" H 8000 1700 50  0001 C CNN
	4    7950 1500
	1    0    0    -1  
$EndComp
$EndSCHEMATC
