# FK6 地平星图

这是一个以 FK6 恒星表为核心的交互式地平星图。程序将 FK6 的 J2000
位置按自行传播至启动时刻；之后每秒刷新观测地点的 Alt/Az 坐标。缩放、
平移和星等上限只重投影已缓存的坐标，不会重复进行自行传播。

## 目录与数据

启动程序需要下列文件：

```text
data/fk6_stars.bin
config/observatory.json
renderer/projection.py
main.py
```

仓库已包含 `data/fk6_stars.bin`。若需要从两份 FK6 原始 ECSV 重建它：

```powershell
E:\Anaconda\envs\tianguangso\python.exe tools\convert_fk6.py
```

`python main.py` 只读取本地 `data/fk6_stars.bin`，不会下载 FK6 原始数据。
太阳、月亮和行星计算使用 Astropy；程序已关闭 Astropy 的 IERS 自动联网刷新，
避免每次启动时尝试下载运行时数据。若手动运行 `tools\download_fk6.py`，
脚本会优先检查 `data\raw` 中已有 ECSV，行数正确时直接复用本地文件。

## 环境

建议使用项目创建的 Anaconda 环境：

```powershell
conda activate tianguangso
```

所需 Python 包：`numpy`、`pygame`、`astropy`。若环境尚未安装：

```powershell
pip install numpy pygame astropy
```

## 启动

在 `FK6\python` 目录中运行：

```powershell
python main.py
```

或使用环境的解释器：

```powershell
E:\Anaconda\envs\tianguangso\python.exe main.py
```

默认观测地点从 [`config/observatory.json`](config/observatory.json) 读取，包含东经为正的
`longitude_deg`、北纬为正的 `latitude_deg` 与 WGS84 椭球高 `height_m`。MFC 启动时会将
界面中的地点值显式传给 Python，因此不再依赖代码内置地点。

### 常用参数

```powershell
# 指定观测地点（覆盖配置文件）
python main.py --longitude 121.55 --latitude 29.87 --height 15

# 指定星等上限
python main.py --magnitude-limit 7.0

# 固定时刻（适合复现实验或导出）
python main.py --time "2026-08-07T20:00:00+08:00"

# 输出一张 PNG 而不打开窗口
python main.py --time "2026-08-07T20:00:00+08:00" --export rendered\sky.png

# 将每次 GOTO/Track 的目标坐标原子写入给宿主程序读取
python main.py --goto-output $env:TEMP\starmap_goto.json
```

使用 `python main.py --help` 可查看全部参数。

## 操作

| 操作 | 功能 |
| --- | --- |
| 鼠标滚轮 | 缩放视野 |
| 鼠标悬停星体 | 显示星体信息；恒星包含 FK6 字段，月球包含月相和亮面比例 |
| 鼠标左键点击 Star / Moon / Sun / Planet | 打开统一的 Selected Object 详情面板 |
| Selected Object → `GOTO` | 生成该对象当前 Alt/Az 目标；黄色准星是目标，绿色准星是实际望远镜位置 |
| Selected Moon / Planet → `TRACK` | 每秒更新其 GOTO 目标并向宿主回传；再次点击 `STOP TRACK` 停止 |
| 拖动窗口边缘 / 角落 | 调整 Pygame 绘图窗口大小并自动重排星图 |
| `W/A/S/D` 或方向键 | 平移完整星图视图；地平网格、方位标签与天球图层同步移动 |
| `[` / `]` | 减少 / 增加 Vmag 上限 |
| `0` | 重置缩放与平移 |
| `1` | 开关主要星座 |
| `2` | 开关赤道坐标网格、赤纬线、天赤道、黄道和 NCP |
| `3` | 开关 Alt/Az 地平坐标辅助网格和 Zenith；最外围地平线圆仍保留 |
| `4` | 开关完整 24h Hour Angle 时角线 |
| `Esc` | 退出 |

## 图层与数据来源

- FK6：4,150 颗恒星，保存 RA、Dec、pmRA*、pmDE、Vmag。
- 主要星座：内置 J2000 坐标线，不依赖 Hipparcos 文件。
- 太阳、月亮和行星：使用 Astropy 内置历表，地平线以下对象不绘制。
  月球从普通对象层拆出，单独绘制为带亮面比例和月相百分比的小圆盘。
- 恒星显示：Vmag 仅决定星点亮度和大小；星等上限只决定是否显示。
  鼠标悬停恒星时，信息框基于当前屏幕坐标命中检测；`FK6 index`
  是合并二进制文件中的记录序号。
- 选择模型：`SelectedObject` 是 Star、Moon、Sun 和 Planet 的共用选择模型；新天体
  只需创建该对象，而不需再加一套点击与 GOTO 状态。
- 望远镜位置：`actual_telescope_position` 与 `goto_target` 是独立状态。MFC 的
  Actual Altitude/Azimuth 只定义绿色实际位置；GOTO 不会伪装成设备已经移动。
- MFC 回传：MFC 为每次启动传入独立的 `--goto-output` 临时 JSON 文件。Python 用
  临时文件替换的方式原子发布 `sequence`、对象和 Alt/Az；MFC 每 250 ms 读取新序列并显示
  在独立的 **GOTO target returned from Python** 区域，适合日后替换为真实设备命令。
- 图层开关：恒星必须显示；星座、赤道坐标网格/黄道、Alt/Az 网格和 Hour Angle
  时角线是可开关的辅助层。
- 坐标辅助线：Alt/Az 网格包含 15°、30°、45°、60°、75° 高度圈；
  Hour Angle 网格覆盖完整 24h。
- 赤道坐标网格包含 Dec -60°、-30°、0°、+30°、+60° 赤纬线、
  NCP 和 Celestial Equator；其中 Dec 0° 即天赤道。
- 地平线圆和 N/E/S/W 是星图基础框架，关闭 Alt/Az 网格时仍会绘制；
  Zenith 跟随 Alt/Az 网格一起开关。

绘制顺序为：地平网格、时角线、天赤道/黄道/星座、恒星、北天极和太阳系天体、HUD。

缩放与平移是统一的视图变换，因此会同时作用于地平线圆、坐标网格、
N/E/S/W、天顶及全部天球图层；时间刷新只更新天文坐标，不改变当前视图。

## 测试

```powershell
python -m unittest discover -s tests -p "test_fk6_*.py" -v
```

测试覆盖 FK6 二进制读取、批量自行传播、Alt/Az 转换、坐标缓存、平移边界、
辅助图层开关、文字防重叠、统一对象模型，以及 Moon 的点击 → GOTO → JSON 坐标回传和
Track Object 更新流程。
