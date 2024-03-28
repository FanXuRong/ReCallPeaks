<a name="Ujy60"></a>
# ReCallPeaks
ReCallPeaks is a Python project for peak detection in sequencing data, relying on samtools and bamCoverage tools for data manipulation.
<a name="Dependencies"></a>
## Dependencies
The execution of this project demands the support of the following libraries or modules:

1. pandas
2. glob2
3. re
4. numpy
5. multiprocessing
6. os
7. scipy
8. math
9. optparse

Plus, the tools to be installed in the environment are:

- samtools
- bamCoverage

Please make sure the mentioned dependencies are installed in your environment before starting out.
<a name="9b2a4b69"></a>
## Usage Guide

1.  Clone the ReCallPeaks code repository 
```bash
git clone https://github.com/XuechengFang/ReCallPeaks.git
cd ReCallPeaks
```

2.  (Optional) Create and activate a Python virtual environment 
```bash
python3 -m venv env
source env/bin/activate
```

3.  Install Python dependencies 
```bash
pip install pandas glob2 scipy
```

4.  Set the running environment with samtools and bamCoverage in the PATH. 
5.  Prepare input data, .bedgraph and .bam files. 
6.  Run the `ReCallPeaks.py` script via Python or the command line interface. This script accepts the following parameters: 
   - `--bedgraphDirectory`: The directory path where the Bedgraph files reside.
   - `--binsize_length`: (Optional) Slicing window length, default is 100.
   - `--threads`: (Optional) The number of threads used by the program, default is 1.
   - `--control_prefix`: The prefix of a control file.
   - `--outRessultcsv`: (Optional) The name of the output file, default is "ReCallPeaks_Results.csv".
   - `--peaks_height_min`: (Optional) The minimum peak height when peak hunting, default is 20.
   - `--peaks_height_max`: (Optional) The maximum peak height when peak hunting, default is 1000.
   - `--peaks_spacing`: (Optional) The peak spacing value when peak hunting, default is 100.
   - `--threshold`: (Optional) The threshold for the multiple differences between the experimental and control groups, default is 2.
   - `--C_peaks_length`: (Optional) The length of the candidate peak area, default is 1000.

Running example: 
```bash
python ReCallPeaks.py --bedgraphDirectory . --binsize_length 100 --threads 5 --control_prefix TB0135 --peaks_height_min 20 --peaks_height_max 700 --peaks_spacing 100 --threshold 2 --C_peaks_length 1000
```

---

<a name="ReCallPeaks"></a>
# ReCallPeaks
ReCallPeaks是一个用于测序数据进行峰值检测的Python项目，依赖于samtools和bamCoverage工具进行数据处理。
<a name="6860b943"></a>
## 依赖
该项目的运行需要以下工具库或模块的支持：

1. pandas
2. glob2
3. re
4. numpy
5. multiprocessing
6. os
7. scipy
8. math
9. optparse

以及已经安装在环境中的工具：

- samtools
- bamCoverage

在开始之前，请确保您的环境中已安装了这些依赖。
<a name="981f67ed"></a>
## 使用指南

1.  克隆ReCallPeaks代码库 
```bash
git clone https://github.com/XuechengFang/ReCallPeaks.git
cd ReCallPeaks
```
 

2.  （可选）创建并激活Python虚拟环境 
```bash
python3 -m venv env
source env/bin/activate
```
 

3.  安装Python依赖 
```bash
pip install pandas glob2 scipy
```
 

4.  设置运行环境，确保samtools和bamCoverage在PATH路径中。 
5.  准备输入数据，.bedgraph和.bam文件。 
6.  通过Python或命令行界面运行`ReCallPeaks.py`脚本。此脚本接受以下参数： 
   - `--bedgraphDirectory`：Bedgraph文件所在的目录路径。
   - `--binsize_length`： （可选）切片窗口长度，默认为100。
   - `--threads`： （可选）程序使用的线程数量，默认为1。
   - `--control_prefix`：控制组文件的前缀。
   - `--outRessultcsv`： （可选）输出文件名，默认为"ReCallPeaks_Results.csv"。
   - `--peaks_height_min`： （可选）寻峰时的最小高度，默认为20。
   - `--peaks_height_max`： （可选）寻峰时的最大高度，默认为1000。
   - `--peaks_spacing`： （可选）寻峰时的间隔值，默认为100。
   - `--threshold`： （可选）实验组与对照组的倍数差异阈值，默认为2。
   - `--C_peaks_length`： （可选）候选峰值区域长度，默认为1000。

运行示例： 
```bash
python ReCallPeaks.py --bedgraphDirectory . --binsize_length 100 --threads 5 --control_prefix TB0135 --peaks_height_min 20 --peaks_height_max 700 --peaks_spacing 100 --threshold 2 --C_peaks_length 1000
```
查看所有可用选项，可以运行 `python ReCallPeaks.py --help`。 
