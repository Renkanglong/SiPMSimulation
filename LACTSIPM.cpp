#include "LACTSIPM.h"
#include <iostream>


std::string LACTSIPM::WAVEFROMSAMPLEFILE_PATH = "file/Waveform_HG.csv";
std::string LACTSIPM::BASELINESAMPLEFILE_PATH = "file/Baseline_HG.csv";
std::string LACTSIPM::PDEDISFILE_PATH = "file/PDE_ARRAY.csv";
std::string LACTSIPM::CHERENKOVDISFILE_PATH = "file/Clw.csv";
std::vector<double> LACTSIPM::WAVEFORM_TIME;
std::vector<double> LACTSIPM::WAVEFORM_COEFFICIENT;
std::vector<double> LACTSIPM::CONVOLUTION_TIME;
std::vector<double> LACTSIPM::CONVOLUTION_COEFFICIENT;
std::vector<double> LACTSIPM::PDE_ARRAY;
std::vector<double> LACTSIPM::CHERENKOV_WAVELENGTH_ARRAY;
std::vector<double> LACTSIPM::CHERENKOV_FREQUENCE_ARRAY;
bool LACTSIPM::NSB_ON = true;
bool LACTSIPM::DCR_ON = true;
bool LACTSIPM::PDE_ON = true;
bool LACTSIPM::CROSSTALK_ON = true;
bool LACTSIPM::AFTERPULSE_ON = true;
bool LACTSIPM::ELENOISE_ON = true;
int LACTSIPM::SPAD_NUM_X = 600;
int LACTSIPM::SPAD_NUM_Y = 600;
int LACTSIPM::TIME_WIN = 256;
int LACTSIPM::CONVOLUTION_WIN = 256;
std::random_device LACTSIPM::SEED;
std::string LACTSIPM::GAIN_TYPE = "HG";
std::string LACTSIPM::CROSSTALK_MODEL = "Binomial";
std::string LACTSIPM::WAVEFORM_UNIT = "mV";
int LACTSIPM::SAMPLE_INTERVAL = 1;
double LACTSIPM::NORM_THRESHOLD = 0.05;
std::mt19937 LACTSIPM::PROBABILITY_RANDOM1{ SEED() }; 
std::mt19937 LACTSIPM::PROBABILITY_RANDOM2{ SEED() }; 
std::mt19937 LACTSIPM::PROBABILITY_RANDOM3{ SEED() }; 
double LACTSIPM::CROSSTALK_PROBABILITY = 0.0173;
double LACTSIPM::AFTERPULSE_PROBABILITY = 0.01;
double LACTSIPM::AFTERPULSE_TEXP = 20;
double LACTSIPM::AFTERPULSE_TREC = 16.7;
double LACTSIPM::DCR_RATE = 8.448e-3;
double LACTSIPM::NSB_RATE = 0.12;
double LACTSIPM::PE_MV_RATIO = 3.09;


LACTSIPM::LACTSIPM()
{
}

LACTSIPM::LACTSIPM(const int& x, const int& y)
{
	SPAD_NUM_X = x;
	SPAD_NUM_Y = y;
}

LACTSIPM::~LACTSIPM()
{
}

template<typename T>
void LACTSIPM::setParameter(const std::string& parametername, T value)
{
	if (parametername == "NSB") NSB_ON = static_cast<bool>(value);
	else if (parametername == "DCR") DCR_ON = static_cast<bool>(value);
	else if (parametername == "PDEResponse") PDE_ON = static_cast<bool>(value);
	else if (parametername == "CrossTalk") CROSSTALK_ON = static_cast<bool>(value);
	else if (parametername == "AfterPulse") AFTERPULSE_ON = static_cast<bool>(value);
	else if (parametername == "ElecNoise") ELENOISE_ON = static_cast<bool>(value);
	else if (parametername == "CrosstalkPro") CROSSTALK_PROBABILITY = static_cast<double>(value);
	else if (parametername == "AfterpulsePro") AFTERPULSE_PROBABILITY = static_cast<double>(value);
	else if (parametername == "AfterpulseTexp") AFTERPULSE_TEXP = static_cast<double>(value);
	else if (parametername == "AfterpulseTrec") AFTERPULSE_TREC = static_cast<double>(value);
	else if (parametername == "DCRRate") DCR_RATE = static_cast<double>(value);
	else if (parametername == "CrosstalkModel") CROSSTALK_MODEL = value;
	else std::cerr << "Undefined Parameter:" << parametername << "\n";
}


void LACTSIPM::initEvent(const std::vector<std::pair<double, double>>& pos, const std::vector<double>& arrive_time, const std::vector<double>& pe)
{
	const double cell_sizex = 24.4 / SPAD_NUM_X;
	const double cell_sizey = 24.4 / SPAD_NUM_Y;
	for (int i = 0;i < pos.size();++i)
	{
		double x = pos[i].first;
		double y = pos[i].second;
		int ix = static_cast<int>(x / cell_sizex);
		int iy = static_cast<int>(y / cell_sizey);
		Event.pos[i] = {ix,iy};
	}
	Event.t_pos = arrive_time;
	Event.pe = pe;
	Event.nfired = pos.size();
	Event.pos_map.resize(SPAD_NUM_X, std::vector<bool>(SPAD_NUM_Y, false));
	for (auto& pos : Event.pos)
	{
		Event.pos_map[pos.first][pos.second] = true;
	}
}

void LACTSIPM::initEvent()
{
	Event.nfired = 0;
	Event.pos_map.resize(SPAD_NUM_X, std::vector<bool>(SPAD_NUM_Y, false));
}

template<typename T>
void LACTSIPM::setWaveformFormat(const std::string& parametername, T value)
{
	if (parametername == "GainType") GAIN_TYPE = value;
	else if (parametername == "WaveformUnit") WAVEFORM_UNIT = value;
	else if (parametername == "SampleInterwal") SAMPLE_INTERVAL = static_cast<int>(value);
	else if (parametername == "TimeWin") TIME_WIN = static_cast<int>(value);
	else if (parametername == "CovolutionWin") CONVOLUTION_WIN = static_cast<int>(value);
	else if (parametername == "NormThreshold") NORM_THRESHOLD = static_cast<double>(value);
	else std::cerr << "Undefined parameter: " << parametername << "\n";
}

void LACTSIPM::setSamplingWaveformPath(const std::string& filepath)
{
	if (GAIN_TYPE == "HG") WAVEFROMSAMPLEFILE_PATH = filepath + "/Waveform_HG.csv";
	if (GAIN_TYPE == "LG") WAVEFROMSAMPLEFILE_PATH = filepath + "/Waveform_LG.csv";
	if (GAIN_TYPE == "HG") BASELINESAMPLEFILE_PATH = filepath + "/Baseline_HG.csv";
	if (GAIN_TYPE == "LG") BASELINESAMPLEFILE_PATH = filepath + "/Baseline_LG.csv";
	PDEDISFILE_PATH = filepath + "/PDE_ARRAY.csv";
	CHERENKOVDISFILE_PATH = filepath + "/Clw.csv";
}

void LACTSIPM::generateCrossTalk()
{
	uint8_t crossstate = 0;
	int n = 4;
	int max_cascade = 5;
	int row, col, ctk_idx, ctk_count;
	int ctk_estimate = Event.nfired * n * (max_cascade + 1);
	std::vector<std::pair<int, int>> new_trigger, ctk_trigger;
	new_trigger.assign(Event.pos.begin(), Event.pos.begin() + Event.nfired);
	Event.pos.resize(ctk_estimate);
	Event.t_pos.resize(ctk_estimate);
	Event.pe.resize(ctk_estimate);
	int total_added = 0; 
	for (int cascade = 0;cascade <= max_cascade;++cascade)
	{
		ctk_idx = Event.nfired + total_added;
		ctk_count = 0;
		ctk_trigger.clear();
		ctk_trigger.resize(new_trigger.size() * n);
		for (int i = 0;i < new_trigger.size();++i)
		{
			int ctk_n = crossNum(n, CROSSTALK_MODEL);
			if (ctk_n <= 0) continue;
			crossstate = crossState(ctk_n, n);
			row = new_trigger[i].first;
			col = new_trigger[i].second;
			if (((crossstate & (1 << 0)) != 0 && col - 1 >= 0))
			{
				Event.pos_map[row][col - 1] = true;
				Event.pos[ctk_idx] = { row, col - 1 };
				Event.t_pos[ctk_idx] = Event.t_pos[i];
				Event.pe[ctk_idx] = 1.0;
				ctk_trigger[ctk_count++] = { row, col - 1 };
				++ctk_idx;
			}
			if (((crossstate & (1 << 1)) != 0 && row - 1 >= 0))
			{
				Event.pos_map[row - 1][col] = true;
				Event.pos[ctk_idx] = { row - 1, col };
				Event.t_pos[ctk_idx] = Event.t_pos[i];
				Event.pe[ctk_idx] = 1.0;
				ctk_trigger[ctk_count++] = { row - 1, col };
				++ctk_idx;
			}
			if (((crossstate & (1 << 2)) != 0 && col + 1 < Event.pos_map[0].size() - 1))
			{
				Event.pos_map[row][col + 1] = true;
				Event.pos[ctk_idx] = { row , col + 1 };
				Event.t_pos[ctk_idx] = Event.t_pos[i];
				Event.pe[ctk_idx] = 1.0;
				ctk_trigger[ctk_count++] = { row, col + 1 };
				++ctk_idx;
			}
			if (((crossstate & (1 << 3)) != 0 && row + 1 <= Event.pos_map.size() - 1))
			{
				Event.pos_map[row + 1][col] = true;
				Event.pos[ctk_idx] = { row + 1, col };
				Event.t_pos[ctk_idx] = Event.t_pos[i];
				Event.pe[ctk_idx] = 1.0;
				ctk_trigger[ctk_count++] = { row + 1, col };
				++ctk_idx;
			}
		}
		if (ctk_count == 0) break;
		ctk_trigger.resize(ctk_count);
		new_trigger = std::move(ctk_trigger);
		total_added += ctk_count;
	}
	Event.nfired += total_added;
	Event.pos.resize(Event.nfired);
	Event.t_pos.resize(Event.nfired);
	Event.pe.resize(Event.nfired);
}

void LACTSIPM::generateDCR()
{
	std::poisson_distribution<> dis(DCR_RATE * TIME_WIN);
	std::uniform_int_distribution<> dis1(.0, SPAD_NUM_X - 1);
	std::uniform_int_distribution<> dis2(.0, SPAD_NUM_Y - 1);
	std::uniform_real_distribution<> dis3(0.0, TIME_WIN);
	int dcr_num = dis(PROBABILITY_RANDOM1);
	int dcr_row, dcr_col;
	double dcr_time;
	Event.pos.resize(Event.nfired + dcr_num);
	Event.t_pos.resize(Event.nfired + dcr_num);
	Event.pe.resize(Event.nfired + dcr_num);
	for (int i = 0;i < dcr_num;++i)
	{
		dcr_row = dis1(PROBABILITY_RANDOM2);
		dcr_col = dis2(PROBABILITY_RANDOM2);
		dcr_time = dis3(PROBABILITY_RANDOM3);
		Event.pos[Event.nfired + i] = { dcr_row, dcr_col };
		Event.t_pos[Event.nfired + i] = dcr_time;
		Event.pe[Event.nfired + i] = 1.0;
		Event.pos_map[dcr_row][dcr_col] = true;
	}
	Event.nfired = Event.nfired + dcr_num;
}

void LACTSIPM::generateNSB()
{
	std::poisson_distribution<> dis(DCR_RATE * TIME_WIN);
	std::uniform_int_distribution<> dis1(0, SPAD_NUM_X - 1);
	std::uniform_int_distribution<> dis2(0, SPAD_NUM_Y - 1);
	std::uniform_real_distribution<> dis3(0.0, TIME_WIN);
	int nsb_num = dis(PROBABILITY_RANDOM1);
	int nsb_row, nsb_col;
	double nsb_time;
	Event.pos.resize(Event.nfired + nsb_num);
	Event.t_pos.resize(Event.nfired + nsb_num);
	Event.pe.resize(Event.nfired + nsb_num);
	for (int i = 0;i < nsb_num;++i)
	{
		nsb_row = dis1(PROBABILITY_RANDOM2);
		nsb_col = dis2(PROBABILITY_RANDOM2);
		nsb_time = dis3(PROBABILITY_RANDOM3);
		Event.pos[Event.nfired + i] = { nsb_row, nsb_col };
		Event.t_pos[Event.nfired + i] = nsb_time;
		Event.pe[Event.nfired + i] = 1.0;
		Event.pos_map[nsb_row][nsb_col] = true;
	}
	Event.nfired = Event.nfired + nsb_num;
}

void LACTSIPM::generateElecNoise()
{
	static bool has_computed = 0;
	static double mean = 0;
	static double stdev = 0;
	static std::vector<double> data;

	if (!has_computed)
	{
		rapidcsv::Document doc(BASELINESAMPLEFILE_PATH, rapidcsv::LabelParams(0, -1));
		data.reserve(doc.GetColumnCount() * doc.GetRowCount());
		for (int i = 0;i < doc.GetColumnCount();++i)
		{
			const std::vector<double>& event = doc.GetColumn<double>(i);
			data.insert(data.end(), event.begin(), event.end());
		}
		stdev = TMath::RMS(data.begin(), data.end());
		mean = TMath::Mean(data.begin(), data.end());
		has_computed = true;
	}
	std::normal_distribution<> Gaussian(mean, stdev);
	for (int i = 0;i < CONVOLUTION_COEFFICIENT.size();++i)
	{
		CONVOLUTION_COEFFICIENT[i] += Gaussian(SEED);
	}
}

void LACTSIPM::Convolution()
{
	readSamplingWfmData();
	const double waveform_interval = WAVEFORM_TIME[1] - WAVEFORM_TIME[0];
	const size_t n_bins = static_cast<size_t>(CONVOLUTION_WIN / waveform_interval);
	CONVOLUTION_TIME.resize(n_bins, 0.0);
	CONVOLUTION_COEFFICIENT.assign(n_bins, 0.0);
	for (int i = 0; i < n_bins; ++i)
		CONVOLUTION_TIME[i] = (i + 1) * waveform_interval;
	std::vector<double> pe_sum(n_bins, 0.0);
	for (int i = 0; i < Event.nfired; ++i)
	{
		double t = Event.t_pos[i];
		double pe = Event.pe[i];
		if (t < 1.0 || t >= CONVOLUTION_WIN) continue;
		size_t bin_idx = static_cast<size_t>(std::floor((t) / waveform_interval)) - 1;
		if (bin_idx >= n_bins) continue;
		pe_sum[bin_idx] += Event.pe[i];
	}
	for (int i = 0; i < n_bins; ++i)
	{
		if (pe_sum[i] == 0.0) continue;
		for (int j = 0; j < WAVEFORM_TIME.size(); ++j)
		{
			if (i + j >= n_bins) continue;
			CONVOLUTION_COEFFICIENT[i + j] += WAVEFORM_COEFFICIENT[j] * pe_sum[i] * PE_MV_RATIO;
		}
	}
}

void LACTSIPM::readPDEData()
{
	static bool has_read_pde = 0;
	static bool has_read_Cherenkov_wavelength = 0;
	if (!has_read_pde)
	{
		rapidcsv::Document doc(PDEDISFILE_PATH, rapidcsv::LabelParams(0, -1));
		PDE_ARRAY = doc.GetColumn<double>(1);
		has_read_pde = true;
	}
	if (!has_read_Cherenkov_wavelength)
	{
		rapidcsv::Document doc1(CHERENKOVDISFILE_PATH, rapidcsv::LabelParams(-1, -1));
		CHERENKOV_WAVELENGTH_ARRAY = doc1.GetColumn<double>(0);
		CHERENKOV_FREQUENCE_ARRAY = doc1.GetColumn<double>(1);
		has_read_Cherenkov_wavelength = true;
	}
}

void LACTSIPM::readSamplingWfmData()
{
	static bool has_read_wave = 0;
	if (!has_read_wave)
	{
		rapidcsv::Document doc(WAVEFROMSAMPLEFILE_PATH, rapidcsv::LabelParams(0, -1));
		const std::vector<double>& waveform_coefficient_samlple = doc.GetColumn<double>(doc.GetColumnCount() - 1);
		int point = waveform_coefficient_samlple.size() / SAMPLE_INTERVAL;
		WAVEFORM_COEFFICIENT.resize(point + 1);
		int idx = 0;
		for (int i = 0;i < waveform_coefficient_samlple.size();i += SAMPLE_INTERVAL)
		{
			WAVEFORM_COEFFICIENT[idx++] = waveform_coefficient_samlple[i];
		}
		waveformNorm();
		WAVEFORM_TIME.resize(WAVEFORM_COEFFICIENT.size());
		for (int i = 0;i < WAVEFORM_COEFFICIENT.size();++i)
		{
			WAVEFORM_TIME[i] = i;
		}
		has_read_wave = true;
	}
}

void LACTSIPM::generateAfterPulse()
{
	std::uniform_real_distribution<> dis(0.0, 1.0);
	std::exponential_distribution<double> dis1(1.0 / AFTERPULSE_TEXP);
	int aft_row, aft_col;
	double p_random, aft_time, aft_pe, exp_time;
	Event.pos.resize(Event.nfired * 2);
	Event.t_pos.resize(Event.nfired * 2);
	Event.pe.resize(Event.nfired * 2);
	int aft_idx = Event.nfired;
	for (int i = 0;i < Event.nfired;++i)
	{
		aft_row = Event.pos[i].first;
		aft_col = Event.pos[i].second;
		p_random = dis(PROBABILITY_RANDOM1);
		if (p_random < AFTERPULSE_PROBABILITY)
		{
			exp_time = 1 + dis1(PROBABILITY_RANDOM2);
			aft_time = Event.t_pos[i] + exp_time;
			if (exp_time < TIME_WIN)
			{
				aft_pe = 1 - std::exp(-exp_time/AFTERPULSE_TREC);
				Event.pos[aft_idx] = { aft_row, aft_col };
				Event.t_pos[aft_idx] = aft_time;
				Event.pe[aft_idx] = aft_pe;
				Event.pos_map[aft_row][aft_col] = true;
				++aft_idx;
			}
		}
	}
	Event.nfired = aft_idx;
	Event.pos.resize(Event.nfired);
	Event.t_pos.resize(Event.nfired);
	Event.pe.resize(Event.nfired);
}

std::vector<double> LACTSIPM::getWaveform_time()
{
	return WAVEFORM_TIME;
}

std::vector<double> LACTSIPM::getWaveform_coefficient()
{
	return WAVEFORM_COEFFICIENT;
}

std::vector<double> LACTSIPM::getConvolution_time()
{
	return CONVOLUTION_TIME;
}

std::vector<double> LACTSIPM::getConvolution_coefficient()
{
	return CONVOLUTION_COEFFICIENT;
}

double LACTSIPM::generateSignal(const std::vector<std::pair<double, double>>& pos, const std::vector<double>& arrive_time, const std::vector<double>& pe)
{
	initEvent(pos, arrive_time, pe);
	if (NSB_ON) generateNSB();
	if (PDE_ON) generatePDEResponse();
	if (DCR_ON) generateDCR();
	if (CROSSTALK_ON) generateCrossTalk();
	if (AFTERPULSE_ON) generateAfterPulse();
	reorder();
	Convolution();
	if (ELENOISE_ON) generateElecNoise();
	setWaveformUnit();
	double Cout = integratedCharge();
	return Cout;
}


inline int LACTSIPM::Combination(int n, int k)
{
	if (k > n) return 0;
	if (k == n || k == 0) return 1;
	int val = 1;
	for (int i = 1; i <= k; ++i)
	{
		val *= (n - i + 1);
		val /= i;
	}
	return val;
}

inline int LACTSIPM::crossNum(int n, std::string type)
{
	if (type == "Binomial")
	{
		double p = 1 - pow(1 - CROSSTALK_PROBABILITY, 0.25);
		std::uniform_real_distribution<> dis(.0, 1.0);
		double p_random = dis(PROBABILITY_RANDOM1);
		double cross_binomial = 0;
		for (int i = n ; i >= 1; --i)
		{
			cross_binomial = Combination(n, i) * pow(p, i) * pow(1 - p, n - i);
			if (p_random < cross_binomial) return i;
		}
		return 0;
	}
	if (type == "Poisson")
	{
		std::uniform_real_distribution<> dis(.0, 1.0);
		double p_random = dis(PROBABILITY_RANDOM1);
		double lamda = -log(1 - CROSSTALK_PROBABILITY);
		double cross_possion = 0;
		for (int i = n; i >= 1; --i)
		{
			cross_possion = pow(lamda, i) * exp(-lamda) / tgamma(i + 1);
			if (p_random < cross_possion) return i;
		}
		return 0;
	}
	if (type == "Geometric")
	{
		std::uniform_real_distribution<> dis(.0, 1.0);
		double p_random = dis(PROBABILITY_RANDOM1);
		if (p_random < CROSSTALK_PROBABILITY) return 1;
	}
	return 0;
}

inline uint8_t LACTSIPM::crossState(int n, int State)
{
	uint8_t crossstate = 0;
	std::uniform_int_distribution<> dis(0, State - 1);
	int count = 0;
	while (count < n)
	{
		int bit_pos = dis(PROBABILITY_RANDOM1);
		if ((crossstate & (1 << bit_pos)) == 0)
		{
			crossstate |= (1 << bit_pos); 
			count++;
		}
	}
	return crossstate;
}

double LACTSIPM::max(std::vector<double>& element)
{
	double max_value = element[0];
	for (int i = 1;i < element.size();++i)
	{
		if (element[i] > max_value) max_value = element[i];
	}
	return max_value;
}

void LACTSIPM::generatePDEResponse()
{
	readPDEData();
	int rps_row, rps_col;
	double p_random, p_pde;
	int rps_idx = 0;
	std::uniform_real_distribution<> dis1(.0, 1.0);
	TH1D* Chrenkov_wavelemgth = new TH1D("", "", CHERENKOV_WAVELENGTH_ARRAY.size(), CHERENKOV_WAVELENGTH_ARRAY.front(), CHERENKOV_WAVELENGTH_ARRAY.back());
	for (int i = 0;i < CHERENKOV_WAVELENGTH_ARRAY.size(); ++i) Chrenkov_wavelemgth->SetBinContent(i + 1, CHERENKOV_FREQUENCE_ARRAY[i]);
	for (int i = 0;i < Event.nfired;i++)
	{
		rps_row = Event.pos[i].first;
		rps_col = Event.pos[i].second;
		int wavelength_random = Chrenkov_wavelemgth->GetRandom();
		p_pde = PDE_ARRAY[wavelength_random - 280];
		p_random = dis1(PROBABILITY_RANDOM1);
		if (p_random < p_pde)
		{
			Event.pos[rps_idx] = { rps_row, rps_col };
			Event.t_pos[rps_idx] = Event.t_pos[i];
			Event.pe[rps_idx] = Event.pe[i];
			Event.pos_map[rps_row][rps_col] = true;
			++rps_idx;
		}
	}
	Event.nfired = rps_idx;
	Event.pos.resize(Event.nfired);
	Event.t_pos.resize(Event.nfired);
	Event.pe.resize(Event.nfired);
}

double LACTSIPM::integratedCharge()
{
	const double waveform_interval = WAVEFORM_TIME[1] - WAVEFORM_TIME[0];
	double integratedcharge = 0.0;
	for (double& CONVOLUTION_COEFFICIENT : CONVOLUTION_COEFFICIENT)
		integratedcharge += CONVOLUTION_COEFFICIENT;
	return integratedcharge;
}

inline double LACTSIPM::recoveryAmplitude(double t1, double t2)
{
	double recoveryamplitude = 1 - exp(-(t2 - t1) / AFTERPULSE_TREC);
	return recoveryamplitude;
}

void LACTSIPM::waveformNorm()
{
	double waveform_max = max(WAVEFORM_COEFFICIENT);
	double baseline = 0;
	int leftturningpoint = 0;
	int rightturningpoint = 0;
	int baseline_samples = WAVEFORM_COEFFICIENT.size() / 10;
	baseline = std::accumulate(WAVEFORM_COEFFICIENT.begin(), WAVEFORM_COEFFICIENT.begin() + baseline_samples, 0.0) / baseline_samples;
	for (double& WAVEFORM_COEFFICIENT : WAVEFORM_COEFFICIENT)
		WAVEFORM_COEFFICIENT = (WAVEFORM_COEFFICIENT - baseline) / (waveform_max - baseline);

	for (int i = 0;i < WAVEFORM_COEFFICIENT.size();++i)
	{
		if (WAVEFORM_COEFFICIENT[i] >= NORM_THRESHOLD)
		{
			leftturningpoint = i;
			break;
		}
	}
	for (int i = WAVEFORM_COEFFICIENT.size() - 1;i >= 0;--i)
	{
		if (WAVEFORM_COEFFICIENT[i] >= NORM_THRESHOLD)
		{
			rightturningpoint = i;
			break;
		}
	}
	WAVEFORM_COEFFICIENT.erase(WAVEFORM_COEFFICIENT.begin(), WAVEFORM_COEFFICIENT.begin() + leftturningpoint);
	WAVEFORM_COEFFICIENT.erase(WAVEFORM_COEFFICIENT.begin() + rightturningpoint + 1 - leftturningpoint + 1, WAVEFORM_COEFFICIENT.end());

}

void LACTSIPM::setWaveformUnit()
{
	if (WAVEFORM_UNIT == "Pe")
	{
		for (auto& WAVEFORM_COEFFICIENT : WAVEFORM_COEFFICIENT)
		{
			WAVEFORM_COEFFICIENT /= PE_MV_RATIO;
		}
		for (auto& CONVOLUTION_COEFFICIENT : CONVOLUTION_COEFFICIENT)
		{
			CONVOLUTION_COEFFICIENT /= PE_MV_RATIO;
		}
	}
	if (WAVEFORM_UNIT == "Adc")
	{
		return;
	}
}

void LACTSIPM::setWaveformUnit(const std::string& type)
{
	if (WAVEFORM_UNIT == type)
	{
		return;
	}
	if (WAVEFORM_UNIT == "mV")
	{
		if (type == "Pe")
		{
			for (auto& WAVEFORM_COEFFICIENT : WAVEFORM_COEFFICIENT)
			{
				WAVEFORM_COEFFICIENT /= PE_MV_RATIO;
			}
			for (auto& CONVOLUTION_COEFFICIENT : CONVOLUTION_COEFFICIENT)
			{
				CONVOLUTION_COEFFICIENT /= PE_MV_RATIO;
			}
		}
		if (type == "Adc")
		{
			return;
		}
	}
	if (WAVEFORM_UNIT == "Pe")
	{
		if (type == "mV")
		{
			for (auto& WAVEFORM_COEFFICIENT : WAVEFORM_COEFFICIENT)
			{
				WAVEFORM_COEFFICIENT *= PE_MV_RATIO;
			}
			for (auto& CONVOLUTION_COEFFICIENT : CONVOLUTION_COEFFICIENT)
			{
				CONVOLUTION_COEFFICIENT *= PE_MV_RATIO;
			}
		}
		if (type == "Adc")
		{
			return;
		}
	}
	if (type == "Adc")
	{
		return;
	}
}

void LACTSIPM::reorder()
{
	std::vector<hitEvent> event;
	event.resize(Event.nfired);
	for (int i = 0;i < Event.nfired;++i)
		event[i] = { Event.pos[i],Event.t_pos[i] };
	std::stable_sort(event.begin(), event.end(),
		[](const hitEvent& event1, const hitEvent& event2)
		{
			if (event1.pos != event2.pos)
				return event1.pos < event2.pos;
			return event1.t_pos < event2.t_pos;
		});
	for (int i = 0;i < Event.nfired;++i)
	{
		Event.pos[i] = event[i].pos;
		Event.t_pos[i] = event[i].t_pos;
		if (i > 0 && Event.pos[i] == Event.pos[i - 1])
			Event.pe[i] = recoveryAmplitude(Event.t_pos[i - 1], Event.t_pos[i]);
		else
			Event.pe[i] = 1.0;
	}
}


