#ifndef LACTSIPM_H
#define LACTSIPM_H

#include <vector>
#include <string>
#include <cmath>
#include <random>   // for mt19937
#include <utility>  // for std::pair
#include <cstdint>  // for unit8_t
#include "rapidcsv.h" //for csvread
#include <numeric>  //for accumulate
#include "TMath.h"  //for mean ,rms
#include "TSpline.h" //for CDF
#include "TH1D.h"
#include "TSpectrum.h"
#include "TMultiGraph.h"

class LACTSIPM
{
public:
    LACTSIPM();
    LACTSIPM(const int& x, const int& y);
    ~LACTSIPM();
    struct sipmSTRUCT
    {
        std::vector<std::vector<bool>> pos_map;        // 像素的map（是否触发）
        std::vector<std::pair<int, int>> pos;          // 每个触发像素的xy位置
        std::vector<double> t_pos;                     // 每个触发像素的触发时间
        std::vector<double> pe;                        // 每个触发像素的光电子强度
        int nfired = 0;                                // 总的触发像素数
    };
    struct hitEvent                            
	{                                                  
		std::pair<int, int> pos;
		double t_pos;
        int pe;
	};
    sipmSTRUCT Event;//初始化的一个结构体, 用于记录SiPM的每个Cell上的光子响应和噪声, 后续的串扰、后脉冲和暗噪声都在这一结构体上操作和记录

private:
    static int SPAD_NUM_X;    // SiPM单通道列数
    static int SPAD_NUM_Y;   // SiPM单通道行数
    static std::random_device SEED; //随机种子
    static std::mt19937 PROBABILITY_RANDOM1; // 用于随机概率1
    static std::mt19937 PROBABILITY_RANDOM2; // 用于随机概率2
    static std::mt19937 PROBABILITY_RANDOM3; // 用于随机概率3
    static double CROSSTALK_PROBABILITY;//固定的串扰概率0.0173；
    static double AFTERPULSE_PROBABILITY;//固定的后脉冲概率0.04；
    static double AFTERPULSE_TEXP;//后脉冲指数时间常数20ns；
    static double AFTERPULSE_TREC;//后脉冲恢复时间常数16.7ns;
    static double DCR_RATE;//暗噪声比例,/ns
    static double NSB_RATE;//夜天光比例,/ns
    static std::vector<double> PDE_ARRAY;//存储PDE_ARRAY
    static std::vector<double> CHERENKOV_WAVELENGTH_ARRAY;//存储切伦科夫光波长
    static std::vector<double> CHERENKOV_FREQUENCE_ARRAY;//存储切伦科夫光波长出现频数
    static std::vector<double> WAVEFORM_TIME;  //读取的波形时间
    static std::vector<double> WAVEFORM_COEFFICIENT;//读取的波形系数
    static std::vector<double> CONVOLUTION_TIME;//卷积的时间
    static std::vector<double> CONVOLUTION_COEFFICIENT;//卷积的系数
    static bool NSB_ON;//夜天光开关
    static bool DCR_ON;//暗噪声开关
    static bool PDE_ON;//光电探测效率开关
    static bool CROSSTALK_ON;//光学串扰开关
    static bool AFTERPULSE_ON;//后脉冲开关
    static bool ELENOISE_ON;//电子学噪声开关
    static int TIME_WIN;//产生噪声的时间窗口，256ns
    static int CONVOLUTION_WIN;//卷积输出的时间窗口
    static std::string GAIN_TYPE;//选择高低增益类型，HG、LG
    static std::string CROSSTALK_MODEL;//串扰模型，Geometric、Binomial、Possion
    static std::string WAVEFORM_UNIT;//输出波形单位类型，Pe, mV, Adc
    static std::string WAVEFROMSAMPLEFILE_PATH;   //采样波长文件路径
    static std::string BASELINESAMPLEFILE_PATH;   //采样基线文件路径
    static std::string PDEDISFILE_PATH; //光电探测效率文件路径
    static std::string CHERENKOVDISFILE_PATH; //切伦科夫光波长分布文件路径
    static int SAMPLE_INTERVAL;// 采样间隔，1point/ns
    static double NORM_THRESHOLD;// 归一化阈值，取0.05阈值中间的波形
    static double PE_MV_RATIO; //pe_mv的单位转换关系

public:

    /// <summary>
    /// 产生SiPM的模拟信号. 完成简单模拟时, 仅调用该函数即可按默认SiPM参数和电子学特性产生波形
    /// </summary>
    /// <param name="pos: 光子在SiPM坐标系中的位置坐标"></param>
    /// <param name="arrive_time: 光子到达SiPM的时间"></param>
    /// <param name="pe: 光子数"></param>
    double generateSignal(const std::vector<std::pair<double, double>>& pos, const std::vector<double>& arrive_time, const std::vector<double>& pe);      //产生信号，完成整个模拟过程
    /// <summary>
    /// 噪声参数进行修改
    /// </summary>
    /// <typeparam name=""></typeparam>
    /// <param name="修改的变量名"></param>
    /// <param name="变量值"></param>
    template<typename T>
    void setParameter(const std::string& parmetername, T value);

    /// <summary>
    /// 波形参数进行修改
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="修改的变量名"></param>
    /// <param name="变量值"></param>
    template<typename T>
    void setWaveformFormat(const std::string& parmetername, T value);

    /// <summary>
    /// 初始化SiPM模拟Event结构体
    /// </summary>
    /// <param name="pos"></param>
    /// <param name="arrive_time"></param>
    /// <param name="pe"></param>
    void initEvent(const std::vector<std::pair<double, double>>& pos, const std::vector<double>& arrive_time, const std::vector<double>& pe);
    
    /// <summary>
    /// 用默认参数初始化SiPM模拟Event结构体. 默认无输入光子，初始化为0
    /// </summary>
    void initEvent();
    
    /// <summary>
    /// 生成暗噪声
    /// </summary>
    void generateDCR();
    
    /// <summary>
    /// 生成夜空背景光干扰信号
    /// </summary>
    void generateNSB();
    
    /// <summary>
    /// 根据SiPM的PDE_ARRAY波长特性, 生成SiPM对不同波长光子的响应
    /// </summary>
    void generatePDEResponse();
    
    /// <summary>
    /// 生成光学串扰(不包括延迟串扰)
    /// </summary>
    void generateCrossTalk();   
    
    /// <summary>
    /// 生成后脉冲
    /// </summary>
    void generateAfterPulse();
    
    /// <summary>
    /// 产生电子学噪声
    /// </summary>
    void generateElecNoise();
    
    /// <summary>
    /// 从默认路径读取电子学实测波形
    /// </summary>
    void readSamplingWfmData();    //直接读取成员变量默认文件路径
    
    /// <summary>
    /// 从默认路径读取SiPM的光子探测效率(PDE_ARRAY)-波长曲线
    /// </summary>
    void readPDEData();

    /// <summary>
    /// 设置电子学实测波形文件路径
    /// </summary>
    /// <param name="filepath"></param>
    void setSamplingWaveformPath(const std::string& filepath);


    /// <summary>
    /// 计算卷积的函数
    /// </summary>
    void Convolution();        

    /// <summary>
    /// 计算卷积后的波形面积
    /// </summary>
    double integratedCharge();//卷积波形面积

    /// <summary>
    /// 对相同SiPMPixel光子的时间排序，计算恢复系数
    /// </summary>
    void reorder();

    /// <summary>
    /// 读取变量WAVEFORM_TIME
    /// </summary>
    /// <std::vector<double>></returns>
    std::vector<double> getWaveform_time();

    /// <summary>
    /// 读取变量WAVEFORM_COEFFICIENT
    /// </summary>
    /// <std::vector<double>></returns>
    std::vector<double> getWaveform_coefficient();

    /// <summary>
    /// 读取变量CONVOLUTION_TIME
    /// </summary>
    /// <std::vector<double>></returns>
    std::vector<double> getConvolution_time();

    /// <summary>
    /// 读取变量CONVOLUTION_COEFFICIENT
    /// </summary>
    /// <std::vector<double>></returns>
    std::vector<double> getConvolution_coefficient();

    /// <summary>
    /// 修改输出波形的单位
    /// </summary>
    void setWaveformUnit();

    /// <summary>
    /// 根据parameter，修改输出波形的单位，pe,mv,adc
    /// </summary>
    /// <param name="UNittype"></param>
    void setWaveformUnit(const std::string& type);

private:
    //辅助计算
    
    /// <summary>
    /// 阶乘计算
    /// </summary>
    /// <param name="邻居数(4)"></param>
    /// <param name="k"></param>
    /// <C(n,k)></returns>
    inline int Combination(int n, int k);  

    /// <summary>
    /// 不同分布的串扰个数
    /// </summary>
    /// <param name="邻居数(4)"></param>
    /// <param name="串扰模型(Geometric、Binomial、Possion)）"></param>
    /// <串扰个数></returns>
    inline int crossNum(int n, std::string type);    

    /// <summary>
    /// 位运算给出串扰状态
    /// </summary>
    /// <param name="串扰个数"></param>
    /// <param name="邻居数(4)"></param>
    /// <uint8_t位状态></returns>
    inline uint8_t crossState(int n, int State);  

    /// <summary>
    /// 时间差计算恢复系数
    /// </summary>
    /// <param name="t1"></param>
    /// <param name="t2"></param>
    /// <t1,t2时间差计算的恢复系数></returns>
    inline double recoveryAmplitude(double t1, double t2);

    /// <summary>
    /// 波形归一化
    /// </summary>
    void waveformNorm();

    /// <summary>
    /// 计算vector最大值
    /// </summary>
    /// <param name="element"></param>
    /// <最大值></returns>
    double max(std::vector<double>& element);
};

#endif