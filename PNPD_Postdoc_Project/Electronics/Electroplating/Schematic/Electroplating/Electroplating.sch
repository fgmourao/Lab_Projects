EESchema Schematic File Version 4
EELAYER 29 0
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
L Device:Battery_Cell Bat
U 1 1 5D16A3DD
P 2610 3180
F 0 "Bat" V 2355 3230 50  0000 C CNN
F 1 "9V" V 2446 3230 50  0000 C CNN
F 2 "" V 2610 3240 50  0001 C CNN
F 3 "~" V 2610 3240 50  0001 C CNN
	1    2610 3180
	0    1    1    0   
$EndComp
$Comp
L Connector:Conn_01x02_Female External_Resistor
U 1 1 5D16C19B
P 2110 4180
F 0 "External_Resistor" H 2002 4273 50  0000 C CNN
F 1 "Conn_01x02_Female" H 2002 4274 50  0001 R CNN
F 2 "" H 2110 4180 50  0001 C CNN
F 3 "~" H 2110 4180 50  0001 C CNN
	1    2110 4180
	-1   0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Female External_Ammeter
U 1 1 5D16DE54
P 2110 4680
F 0 "External_Ammeter" H 2002 4773 50  0000 C CNN
F 1 "Conn_01x02_Female" H 2002 4774 50  0001 C CNN
F 2 "" H 2110 4680 50  0001 C CNN
F 3 "~" H 2110 4680 50  0001 C CNN
	1    2110 4680
	-1   0    0    -1  
$EndComp
Wire Wire Line
	2310 4680 2310 4280
$Comp
L Device:LED Led
U 1 1 5D16DFB2
P 2750 4880
F 0 "Led" H 2743 4625 50  0001 C CNN
F 1 "LED" H 2743 4716 50  0000 C CNN
F 2 "" H 2750 4880 50  0001 C CNN
F 3 "~" H 2750 4880 50  0001 C CNN
	1    2750 4880
	-1   0    0    1   
$EndComp
Wire Wire Line
	2310 4780 2310 4880
Wire Wire Line
	2310 4880 2600 4880
$Comp
L Connector:Conn_01x01_Female CH1
U 1 1 5D170C26
P 3510 3780
F 0 "CH1" V 3402 3692 50  0000 R CNN
F 1 "Ch1" V 3357 3692 50  0001 R CNN
F 2 "" H 3510 3780 50  0001 C CNN
F 3 "~" H 3510 3780 50  0001 C CNN
	1    3510 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3510 4480 3510 4080
$Comp
L Connector:Conn_01x01_Female CH2
U 1 1 5D174570
P 3910 3780
F 0 "CH2" V 3802 3692 50  0000 R CNN
F 1 "Ch1" V 3757 3692 50  0001 R CNN
F 2 "" H 3910 3780 50  0001 C CNN
F 3 "~" H 3910 3780 50  0001 C CNN
	1    3910 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3910 4480 3910 3980
$Comp
L Connector:Conn_01x01_Female CH3
U 1 1 5D1750BB
P 4310 3780
F 0 "CH3" V 4202 3692 50  0000 R CNN
F 1 "Ch1" V 4157 3692 50  0001 R CNN
F 2 "" H 4310 3780 50  0001 C CNN
F 3 "~" H 4310 3780 50  0001 C CNN
	1    4310 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	4310 4480 4310 3980
$Comp
L Connector:Conn_01x01_Female CH4
U 1 1 5D179485
P 4710 3780
F 0 "CH4" V 4602 3692 50  0000 R CNN
F 1 "Ch1" V 4557 3692 50  0001 R CNN
F 2 "" H 4710 3780 50  0001 C CNN
F 3 "~" H 4710 3780 50  0001 C CNN
	1    4710 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	4710 4480 4710 3980
$Comp
L Connector:Conn_01x01_Female CH5
U 1 1 5D179492
P 5110 3780
F 0 "CH5" V 5002 3692 50  0000 R CNN
F 1 "Ch1" V 4957 3692 50  0001 R CNN
F 2 "" H 5110 3780 50  0001 C CNN
F 3 "~" H 5110 3780 50  0001 C CNN
	1    5110 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5110 4480 5110 3980
$Comp
L Connector:Conn_01x01_Female CH6
U 1 1 5D17949F
P 5510 3780
F 0 "CH6" V 5402 3692 50  0000 R CNN
F 1 "Ch1" V 5357 3692 50  0001 R CNN
F 2 "" H 5510 3780 50  0001 C CNN
F 3 "~" H 5510 3780 50  0001 C CNN
	1    5510 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5510 4480 5510 3980
$Comp
L Connector:Conn_01x01_Female CH7
U 1 1 5D180C1A
P 5910 3780
F 0 "CH7" V 5802 3692 50  0000 R CNN
F 1 "Ch1" V 5757 3692 50  0001 R CNN
F 2 "" H 5910 3780 50  0001 C CNN
F 3 "~" H 5910 3780 50  0001 C CNN
	1    5910 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5910 4480 5910 3980
$Comp
L Connector:Conn_01x01_Female CH8
U 1 1 5D180C27
P 6310 3780
F 0 "CH8" V 6202 3692 50  0000 R CNN
F 1 "Ch1" V 6157 3692 50  0001 R CNN
F 2 "" H 6310 3780 50  0001 C CNN
F 3 "~" H 6310 3780 50  0001 C CNN
	1    6310 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	6310 4480 6310 3980
$Comp
L Connector:Conn_01x01_Female CH9
U 1 1 5D180C34
P 6710 3780
F 0 "CH9" V 6602 3692 50  0000 R CNN
F 1 "Ch1" V 6557 3692 50  0001 R CNN
F 2 "" H 6710 3780 50  0001 C CNN
F 3 "~" H 6710 3780 50  0001 C CNN
	1    6710 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	6710 4480 6710 3980
$Comp
L Connector:Conn_01x01_Female CH10
U 1 1 5D180C41
P 7110 3780
F 0 "CH10" V 7002 3692 50  0000 R CNN
F 1 "Ch1" V 6957 3692 50  0001 R CNN
F 2 "" H 7110 3780 50  0001 C CNN
F 3 "~" H 7110 3780 50  0001 C CNN
	1    7110 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	7110 4480 7110 3980
$Comp
L Connector:Conn_01x01_Female CH11
U 1 1 5D180C4E
P 7510 3780
F 0 "CH11" V 7402 3692 50  0000 R CNN
F 1 "Ch1" V 7357 3692 50  0001 R CNN
F 2 "" H 7510 3780 50  0001 C CNN
F 3 "~" H 7510 3780 50  0001 C CNN
	1    7510 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	7510 4480 7510 3980
$Comp
L Connector:Conn_01x01_Female CH12
U 1 1 5D180C5B
P 7910 3780
F 0 "CH12" V 7802 3692 50  0000 R CNN
F 1 "Ch1" V 7757 3692 50  0001 R CNN
F 2 "" H 7910 3780 50  0001 C CNN
F 3 "~" H 7910 3780 50  0001 C CNN
	1    7910 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	7910 4480 7910 3980
$Comp
L Connector:Conn_01x01_Female CH13
U 1 1 5D185E85
P 8310 3780
F 0 "CH13" V 8202 3692 50  0000 R CNN
F 1 "Ch1" V 8157 3692 50  0001 R CNN
F 2 "" H 8310 3780 50  0001 C CNN
F 3 "~" H 8310 3780 50  0001 C CNN
	1    8310 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	8310 4480 8310 3980
$Comp
L Connector:Conn_01x01_Female CH14
U 1 1 5D185E92
P 8710 3780
F 0 "CH14" V 8602 3692 50  0000 R CNN
F 1 "Ch1" V 8557 3692 50  0001 R CNN
F 2 "" H 8710 3780 50  0001 C CNN
F 3 "~" H 8710 3780 50  0001 C CNN
	1    8710 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	8710 4480 8710 3980
$Comp
L Connector:Conn_01x01_Female CH15
U 1 1 5D185E9F
P 9110 3780
F 0 "CH15" V 9002 3692 50  0000 R CNN
F 1 "Ch1" V 8957 3692 50  0001 R CNN
F 2 "" H 9110 3780 50  0001 C CNN
F 3 "~" H 9110 3780 50  0001 C CNN
	1    9110 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	9110 4480 9110 3980
$Comp
L Connector:Conn_01x01_Female CH16
U 1 1 5D185EAC
P 9510 3780
F 0 "CH16" V 9402 3692 50  0000 R CNN
F 1 "Ch1" V 9357 3692 50  0001 R CNN
F 2 "" H 9510 3780 50  0001 C CNN
F 3 "~" H 9510 3780 50  0001 C CNN
	1    9510 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	9510 4480 9510 3980
$Comp
L Connector:Conn_01x01_Female REF
U 1 1 5D18736D
P 3510 3380
F 0 "REF" V 3402 3428 50  0000 L CNN
F 1 "Ch1" V 3357 3292 50  0001 R CNN
F 2 "" H 3510 3380 50  0001 C CNN
F 3 "~" H 3510 3380 50  0001 C CNN
	1    3510 3380
	0    1    1    0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C2E
P 6710 4680
F 0 "Push_Botton?" H 6710 4495 50  0001 C CNN
F 1 "Push_Botton" H 6710 4586 50  0000 L CNN
F 2 "" H 6710 4880 50  0001 C CNN
F 3 "" H 6710 4880 50  0001 C CNN
	1    6710 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C21
P 6310 4680
F 0 "Push_Botton?" H 6310 4495 50  0001 C CNN
F 1 "Push_Botton" H 6310 4586 50  0000 L CNN
F 2 "" H 6310 4880 50  0001 C CNN
F 3 "" H 6310 4880 50  0001 C CNN
	1    6310 4680
	0    -1   -1   0   
$EndComp
Connection ~ 6710 4880
Wire Wire Line
	6310 4880 6710 4880
Connection ~ 6310 4880
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D185EA6
P 9510 4680
F 0 "Push_Botton?" H 9510 4495 50  0001 C CNN
F 1 "Push_Botton" H 9510 4586 50  0000 L CNN
F 2 "" H 9510 4880 50  0001 C CNN
F 3 "" H 9510 4880 50  0001 C CNN
	1    9510 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D185E99
P 9110 4680
F 0 "Push_Botton?" H 9110 4495 50  0001 C CNN
F 1 "Push_Botton" H 9110 4586 50  0000 L CNN
F 2 "" H 9110 4880 50  0001 C CNN
F 3 "" H 9110 4880 50  0001 C CNN
	1    9110 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D185E8C
P 8710 4680
F 0 "Push_Botton?" H 8710 4495 50  0001 C CNN
F 1 "Push_Botton" H 8710 4586 50  0000 L CNN
F 2 "" H 8710 4880 50  0001 C CNN
F 3 "" H 8710 4880 50  0001 C CNN
	1    8710 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D185E7F
P 8310 4680
F 0 "Push_Botton?" H 8310 4495 50  0001 C CNN
F 1 "Push_Botton" H 8310 4586 50  0000 L CNN
F 2 "" H 8310 4880 50  0001 C CNN
F 3 "" H 8310 4880 50  0001 C CNN
	1    8310 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C55
P 7910 4680
F 0 "Push_Botton?" H 7910 4495 50  0001 C CNN
F 1 "Push_Botton" H 7910 4586 50  0000 L CNN
F 2 "" H 7910 4880 50  0001 C CNN
F 3 "" H 7910 4880 50  0001 C CNN
	1    7910 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C48
P 7510 4680
F 0 "Push_Botton?" H 7510 4495 50  0001 C CNN
F 1 "Push_Botton" H 7510 4586 50  0000 L CNN
F 2 "" H 7510 4880 50  0001 C CNN
F 3 "" H 7510 4880 50  0001 C CNN
	1    7510 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C3B
P 7110 4680
F 0 "Push_Botton?" H 7110 4495 50  0001 C CNN
F 1 "Push_Botton" H 7110 4586 50  0000 L CNN
F 2 "" H 7110 4880 50  0001 C CNN
F 3 "" H 7110 4880 50  0001 C CNN
	1    7110 4680
	0    -1   -1   0   
$EndComp
Wire Wire Line
	9110 4880 9510 4880
Connection ~ 9110 4880
Connection ~ 8710 4880
Wire Wire Line
	8710 4880 9110 4880
Wire Wire Line
	8310 4880 8710 4880
Connection ~ 8310 4880
Connection ~ 7910 4880
Wire Wire Line
	7910 4880 8310 4880
Wire Wire Line
	7510 4880 7910 4880
Connection ~ 7510 4880
Connection ~ 7110 4880
Wire Wire Line
	7110 4880 7510 4880
Wire Wire Line
	6710 4880 7110 4880
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D180C14
P 5910 4680
F 0 "Push_Botton?" H 5910 4495 50  0001 C CNN
F 1 "Push_Botton" H 5910 4586 50  0000 L CNN
F 2 "" H 5910 4880 50  0001 C CNN
F 3 "" H 5910 4880 50  0001 C CNN
	1    5910 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D179499
P 5510 4680
F 0 "Push_Botton?" H 5510 4495 50  0001 C CNN
F 1 "Push_Botton" H 5510 4586 50  0000 L CNN
F 2 "" H 5510 4880 50  0001 C CNN
F 3 "" H 5510 4880 50  0001 C CNN
	1    5510 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D17948C
P 5110 4680
F 0 "Push_Botton?" H 5110 4495 50  0001 C CNN
F 1 "Push_Botton" H 5110 4586 50  0000 L CNN
F 2 "" H 5110 4880 50  0001 C CNN
F 3 "" H 5110 4880 50  0001 C CNN
	1    5110 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D17947F
P 4710 4680
F 0 "Push_Botton?" H 4710 4495 50  0001 C CNN
F 1 "Push_Botton" H 4710 4586 50  0000 L CNN
F 2 "" H 4710 4880 50  0001 C CNN
F 3 "" H 4710 4880 50  0001 C CNN
	1    4710 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D1750B5
P 4310 4680
F 0 "Push_Botton?" H 4310 4495 50  0001 C CNN
F 1 "Push_Botton" H 4310 4586 50  0000 L CNN
F 2 "" H 4310 4880 50  0001 C CNN
F 3 "" H 4310 4880 50  0001 C CNN
	1    4310 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton?
U 1 1 5D17456A
P 3910 4680
F 0 "Push_Botton?" H 3910 4495 50  0001 C CNN
F 1 "Push_Botton" H 3910 4586 50  0000 L CNN
F 2 "" H 3910 4880 50  0001 C CNN
F 3 "" H 3910 4880 50  0001 C CNN
	1    3910 4680
	0    -1   -1   0   
$EndComp
$Comp
L Switch:SW_Push Push_Botton
U 1 1 5D16EFC7
P 3510 4680
F 0 "Push_Botton" H 3510 4495 50  0001 R CNN
F 1 "Push_Botton" H 3510 4586 50  0000 L CNN
F 2 "" H 3510 4880 50  0001 C CNN
F 3 "" H 3510 4880 50  0001 C CNN
	1    3510 4680
	0    -1   -1   0   
$EndComp
Wire Wire Line
	2900 4880 3510 4880
Wire Wire Line
	3510 4880 3910 4880
Connection ~ 3510 4880
Connection ~ 3910 4880
Wire Wire Line
	3910 4880 4310 4880
Wire Wire Line
	4310 4880 4710 4880
Connection ~ 4310 4880
Connection ~ 4710 4880
Wire Wire Line
	4710 4880 5110 4880
Wire Wire Line
	5110 4880 5510 4880
Connection ~ 5110 4880
Connection ~ 5510 4880
Wire Wire Line
	5510 4880 5910 4880
Wire Wire Line
	5910 4880 6310 4880
Connection ~ 5910 4880
$Comp
L Connector:Conn_01x01_Female CH1_copy
U 1 1 5D1A88DA
P 3110 3780
F 0 "CH1_copy" H 3002 3692 50  0000 C CNN
F 1 "Ch1" V 2957 3692 50  0001 R CNN
F 2 "" H 3110 3780 50  0001 C CNN
F 3 "~" H 3110 3780 50  0001 C CNN
	1    3110 3780
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3110 3980 3110 4080
Wire Wire Line
	3110 4080 3510 4080
Connection ~ 3510 4080
Wire Wire Line
	3510 4080 3510 3980
Wire Wire Line
	2310 3180 2310 4180
$Comp
L Connector:Conn_01x01_Female REF_Copy
U 1 1 5D1AE9C1
P 3110 3380
F 0 "REF_Copy" H 3002 3428 50  0000 L CNN
F 1 "Copy_REF" V 2957 3292 50  0001 R CNN
F 2 "" H 3110 3380 50  0001 C CNN
F 3 "~" H 3110 3380 50  0001 C CNN
	1    3110 3380
	0    1    1    0   
$EndComp
Wire Wire Line
	2610 3180 2810 3180
Wire Wire Line
	3110 3180 3510 3180
Wire Wire Line
	2310 3180 2510 3180
Wire Wire Line
	2810 3180 3110 3180
Connection ~ 2810 3180
Connection ~ 3110 3180
Wire Wire Line
	2510 3180 2510 2780
Connection ~ 2510 3180
Wire Wire Line
	2810 3180 2810 2780
$Comp
L Connector:Conn_01x01_Female External_Bat-
U 1 1 5D1BCBF0
P 2510 2580
F 0 "External_Bat-" H 2402 2492 50  0000 L CNN
F 1 "External_bat-" V 2357 2492 50  0001 R CNN
F 2 "" H 2510 2580 50  0001 C CNN
F 3 "~" H 2510 2580 50  0001 C CNN
	1    2510 2580
	0    -1   -1   0   
$EndComp
$Comp
L Connector:Conn_01x01_Female External_Bat+
U 1 1 5D1BDD49
P 2810 2580
F 0 "External_Bat+" H 2702 2492 50  0000 L CNN
F 1 "Ch1" V 2657 2492 50  0001 R CNN
F 2 "" H 2810 2580 50  0001 C CNN
F 3 "~" H 2810 2580 50  0001 C CNN
	1    2810 2580
	0    -1   -1   0   
$EndComp
$EndSCHEMATC
