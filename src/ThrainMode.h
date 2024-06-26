//**************************************************************************************
//							THRAIN MODE
//  ThrainMode.h
//	 	Defines functionality with finding modes and calculating their properies
//  Reece Boston, Aug, 20, 2023
//**************************************************************************************

#ifndef THRAINMODEH
#define THRAINMODEH

#include <map>
#include <set>
#include <vector>
#include <string>

namespace Calculation {
	struct OutputData;
}

namespace io {
	int write_mode_output(Calculation::OutputData&);
}

namespace mode {

double compare_JCD(double, int, int, double);
double compare_Pekeris(double, int, int, double);

// the following arrays come from tables in JCD-DJM paper:
// 		Christensen-Dalsgaard & Mullan (1994), MNRAS 270, 921-935
// these were very useful for testing in early stages of the code

//columns L = 1, 2, 3 from Table 1 in JCD-DJM paper
const double JCD1_5[3][35] = {
	// L = 1
	{256.8042, 424.9943, 582.9618, 736.2824, 886.9146, 1035.8195, 1183.5437, 1330.4248, 1476.6835, 1622.4704,
		1767.8912, 1913.0229, 2057.9220, 2202.6318, 2347.1851, 2491.6079, 2635.9207, 2780.1399, 2924.2788,
		3068.3483, 3212.3576, 3356.3141, 3500.2242, 3644.0932, 3787.9257, 3931.7257, 4075.4965, 4219.2411,
		4362.9620, 4506.6614, 4650.3412, 4794.0033, 4937.6492, 5081.2801, 5224.8974},
	// L = 2
	{320.2718, 484.2210, 641.7601, 795.5461, 946.8945, 1096.5604, 1245.0155, 1392.5718, 1539.4439, 1685.7838,
		1831.7026, 1977.2825, 2122.5862, 2267.6620, 2412.5474, 2557.2726, 2701.8615, 2846.3336, 2990.7050,
		3131.9888, 3279.1960, 3423.3360, 3567.4165, 3711.4442, 3855.4249, 3999.3634, 4143.2640, 4287.1305,
		4430.9660, 4574.7734, 4718.5552, 4862.3137, 5006.0508, 5149.7683, 5293.4677},
	// L = 3
	{369.0420, 534.0912, 693.1482, 848.4539, 1001.2131, 1152.1531, 1301.7468, 1450.3169, 1598.0915, 1745.2365,
		1891.8754, 2038.1018, 2183.9880, 2329.5905, 2474.9541, 2620.1149, 2765.1019, 2909.9392, 3054.6465,
		3199.2403, 3343.7343, 3488.1403, 3632.4682, 3776.7265, 3920.9226, 4065.0627, 4209.1525, 4353.1996,
		4497.1996, 4641.1649, 4785.0958, 4928.9953, 5072.8660, 5216.7101, 5360.5299}
};

//columns L = 1, 2, 3 from Table 2 in JCD-DJM paper
const double JCD3_0[3][35] = {
	// L = 1
	{337.2152, 463.5718, 590.0694, 716.6289, 843.1066,  969.4331, 1095.5832, 1221.5530, 1347.3484, 1472.9799,
		1598.4594, 1723.7990, 1849.0105, 1974.1051, 2099.0928, 2223.9831, 2348.7844, 2473.5043, 2598.1500,
		2722.7276, 2847.2430, 2971.7011, 3096.1067, 3220.4639, 3344.7767, 3469.0483, 3593.2820, 3717.4807,
		3841.6468, 3965.7827, 4089.8906, 4213.9725, 4338.0301, 4462.0650, 4586.0789},
	// L = 2
	{390.1223, 516.1992, 643.0677, 769.7802, 896.2910, 1022.6135, 1148.7622, 1274.7491, 1400.5849, 1526.2795,
		1651.8424, 1777.2825, 1902.6082, 2027.8277, 2152.9485, 2277.9776, 2402.9216, 2527.7867, 2652.5785,
		2777.3023, 2901.9629, 3026.5648, 3151.1121, 3275.6085, 3400.0577, 3524.4627, 3648.8266, 3773.1522,
		3897.4418, 4021.6979, 4145.9226, 4270.1179, 4394.2856, 4518.4275 ,4642.5450},
	// L = 3
	{428.8391, 558.2981, 686.8732, 814.7027, 942.0267, 1068.9866, 1195.6639, 1322.1090, 1448.3552, 1574.4265,
		1700.3411, 1826.1136, 1951.7561, 2077.2793, 2202.6923, 2328.0032, 2453.2197, 2578.3482, 2703.3949,
		2828.3655, 2953.2651, 3078.0983, 3202.8695, 3327.5828, 3452.2418, 3576.8500, 3701.4106, 3825.9264,
		3950.4002, 4074.8346, 4199.2318, 4323.5941, 4447.9235, 4572.2218, 4696.4908}
};

//columns L = 1, 2, 3 from Table 3 in JCD-DJM paper
const double JCD4_0[3][35] = {
	// L = 1
	{507.0621, 571.8190, 625.9711, 736.0301, 848.1925,  960.7919, 1073.6477, 1186.7036, 1299.9278, 1413.2959,
		1526.7873, 1640.3842, 1754.0710, 1867.8342, 1981.6625, 2095.5458, 2209.4757, 2323.4448, 2437.4470,
		2551.4770, 2665.5302, 2779.6026, 2893.6909, 3007.7923, 3121.9042, 3236.0246, 3350.1515, 3464.2834,
		3578.4189, 3692.5569, 3806.6962, 3920.8360, 4034.9756, 4149.1143, 4263.2516},
	// L = 2
	{648.2038, 711.8195, 791.8287, 874.4382, 933.9077, 1024.6134, 1131.6998, 1242.3945, 1354.4090, 1467.0790,
		1580.1333, 1693.4374, 1806.9155, 1920.5210, 2034.2233, 2148.0012, 2261.8392, 2375.7260, 2489.6527,
		2603.6125, 2717.5999, 2831.6102, 2945.6399, 3059.6858, 3173.7454, 3287.8163, 3401.8968, 3515.9851,
		3630.0799, 3744.1800, 3858.2842, 3972.3917, 4086.5016, 4200.6132, 4314.7259},
	// L = 3
	{713.1345, 782.3745, 833.5703, 936.7684, 996.2500, 1059.2311, 1169.2703, 1282.1336, 1395.5967, 1509.3076,
		1623.1501, 1737.0730, 1851.0500, 1965.0657, 2079.1108, 2193.1790, 2307.2658, 2421.3680, 2535.4830,
		2649.6088, 2763.7438, 2877.8864, 2992.0355, 3106.1901, 3220.3492, 3334.5120, 3448.6777, 3562.8457,
		3677.0153, 3791.1861, 3905.3575, 4019.5291, 4133.7005, 4247.8714, 4362.0413}
};

double freq_JCD(double n, int l, int k);
double omega2_JCD(double n, int l, int k);

// for comparing to known values
double compare_JCD(double n, int l, int k, double w);
double calculate_Pekeris(int l, int k, double Gam1);
double compare_Pekeris(double w, int l, int k, double Gam1);

//this will create the correct mode driver and pass the correct driver type to the mode_finder
int create_modes(Calculation::OutputData &data_out);

template <class MODE>
void save_mode (
	MODE* modetry,
	std::set<int>& kfilled,
	std::map<int,double>& w2filled,
	std::map<int, MODE*>& modefilled
);

void get_min_from_set(const std::set<int>&, const std::map<int,double>&, const int&, int&, double&);
void get_max_from_set(const std::set<int>&, const std::map<int,double>&, const int&, int&, double&);

//this function will find eigenmode and frequency for the desired k,l in the output data
// requires a mode driver type (the <class>)
//  defined below
template <class MODE, class MODEDRIVER> 
int mode_finder(Calculation::OutputData&);

} // namespace mode

//template classes in C++ cannot be split into multiple files, so we must include Mode.pp
#include "ThrainMode.hxx"

#endif