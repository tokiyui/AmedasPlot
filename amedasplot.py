#!/usr/bin/env python
# coding: utf-8

# AMeDAS 観測データプロット図 （地図領域は、緯度・経度のそれぞれ2度の範囲)
# 矢羽,気温(左上),湿球温度（右下）,露点温度(左下),海面更正気圧(右上)
# 黒良先生のプログラムにレーダーGPVを重ね合わせ
# アメダス地点テーブルJSON  https://www.jma.go.jp/bosai/amedas/const/amedastable.json
# アメダス観測データJSON  YYYY/MM/DD HH:mm(JST)  https:https://www.jma.go.jp/bosai/amedas/data/map/{YYYY}{MM}{DD}{HH}{mm}00.json
# 生存圏研究所ダウンロード元サイト  http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original

import argparse, json, math, matplotlib, os, pygrib, pytz, struct, subprocess, sys, csv, re
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import pandas.tseries.offsets as offsets
from datetime import datetime, timedelta
from itertools import repeat
from metpy.units import units
from scipy.interpolate import griddata, interp2d, RectBivariateSpline, RegularGridInterpolator
from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter
from urllib.request import urlopen
import netCDF4 as nc
from ftplib import FTP
import xarray as xr
import requests
from matplotlib.colors import ListedColormap, BoundaryNorm

def dms_to_decimal(dms):
    direction = dms[-1]
    degrees, minutes = dms[:-1].split()
    degrees = float(degrees)
    minutes = float(minutes)
    decimal = degrees + minutes / 60
    if direction in ['S', 'W']:
        decimal *= -1
    return round(decimal,4)
 
def process_url(cou, year, month, day, hour):
    # URLからHTMLを取得
    url = "http://www.meteomanz.com/sy6?cou={}&sh=map4&d1={}&m1={}&y1={}&h1={}Z&l=1".format(cou, day, month, year, hour)
    print(url)
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36'}
    response = requests.get(url, headers= headers)
    html = response.text
 
    # 正規表現パターンを定義
    pattern = r'<span id="(\d+)" style="Display:none;"><table.*?<u>(.*?), (.*?) \((.*?)\/(.*?)\/(.*?)\), (\d+)Z:</u>.*?<td>(.*?)</td>.*?<td>(.*?)</td>.*?<td>(.*?)</td>.*?<td>(.*?)</td>.*?Tmin:(.*?)ºC.*?Tmax:(.*?)ºC.*?</pre>'
 
    # パターンにマッチするデータを抽出
    matches = re.findall(pattern, html, re.DOTALL)
 
    # データをリストに格納
    data_list = []
    for match in matches:
        height = match[5].replace(' m', '')
        temperature = match[7].replace('ÂºC', '')
        if temperature == '':
            temperature = '-'
        humidity = match[8].replace('%', '')
        pressure = '-' if 'N/A' in match[9] or 'm' in match[9] else match[9].replace(' Hpa', '')
       
        if 'calm' in match[10]:
            wind_speed = 'calm'
        else:
            wind_speed = match[10].replace(' Km/h', '') if 'N/A' not in match[10] else '-'
 
        lat = dms_to_decimal(match[3].replace('IN ASIA) (', ''))
        lon = dms_to_decimal(match[4])
 
        wind_info = match[10].split()
        if len(wind_info) == 3:
            wind_direction = wind_info[0]
            wind_speed = wind_info[1]
            wind_speed_decimal = round(float(wind_speed) / 3.6, 1)  # km/hをm/sに変換し、小数点以下1桁までに丸める
        elif wind_info[0] == 'calm':
            wind_direction = 'calm'
            wind_speed_decimal = '-'
        else:
            wind_direction = '-'
            wind_speed_decimal = '-'
 
        data = [match[0], match[1], match[2], lat, lon, height, temperature, humidity, pressure, wind_direction , wind_speed_decimal]
        if lat > 0 and lon > 90 :
          data_list.append(data)
        elif lat > 0 and lon < -150:
          data_list.append(data)
         
    return data_list

## 極大/極小ピーク検出関数                                                             
def detect_peaks(image, filter_size, dist_cut, flag=0):
    # filter_size: この値xこの値 の範囲内の最大値のピークを検出                        
    # dist_cut: この距離内のピークは1つにまとめる                                      
    # flag:  0:maximum検出  0以外:minimum検出                                          
    if flag==0:
      local_max = maximum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant') 
      detected_peaks = np.ma.array(image, mask=~(image == local_max))
    else:
      local_min = minimum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant')
      detected_peaks = np.ma.array(image, mask=~(image == local_min))
    peaks_index = np.where((detected_peaks.mask != True))
    # peak間の距離行例を求める                                                         
    (x,y) = peaks_index
    size=y.size
    dist=np.full((y.size, y.size), -1.0)
    for i in range(size):
      for j in range(size):
        if i == j:
          dist[i][j]=0.0
        elif i>j:
          d = math.sqrt(((y[i] - y[j])*(y[i] - y[j])) + ((x[i] - x[j])*(x[i] - x[j])))
          dist[i][j]= d
          dist[j][i]= d
    # 距離がdist_cut内のpeaksの距離の和と、そのピーク番号を取得する 
    Kinrin=[]
    dSum=[]
    for i in range(size):
      tmpA=[]
      distSum=0.0
      for j in range(size):
        if dist[i][j] < dist_cut and dist[i][j] > 0.0:
          tmpA.append(j)
          distSum=distSum+dist[i][j]
      dSum.append(distSum)
      Kinrin.append(tmpA)
    # Peakから外すPeak番号を求める.  peak間の距離和が最も小さいものを残す              
    cutPoint=[]
    for i in range(size):
      val = dSum[i]
      val_i=image[x[i]][y[i]]
      for k in Kinrin[i]:
        val_k=image[x[k]][y[k]]
        if flag==0 and val_i < val_k:
            cutPoint.append(i)
            break
        if flag!=0 and val_i > val_k:
            cutPoint.append(i)
            break
        if val > dSum[k]:
            cutPoint.append(i)
            break
        if val == dSum[k] and i > k:
            cutPoint.append(i)
            break
    # 戻り値用に外すpeak番号を配列から削除                                             
    newx=[]
    newy=[]
    for i in range(size):
      if (i in cutPoint):
        continue
      newx.append(x[i])
      newy.append(y[i])
    peaks_index=(np.array(newx),np.array(newy))
    return peaks_index

def parse_datetime(arg):
    try:
        # 引数が12桁の数字の場合、YYYYMMDDHHMM形式の文字列を解析
        if len(arg) == 12 and arg.isdigit():
            year = arg[:4]
            month = arg[4:6]
            day = arg[6:8]
            hour = arg[8:10]
            min = arg[10:12]
            dt = datetime(int(year), int(month), int(day), int(hour), int(min))
        elif len(arg) == 10 and arg.isdigit():
            year = arg[:4]
            month = arg[4:6]
            day = arg[6:8]
            hour = arg[8:10]
            min = 0
            dt = datetime(int(year), int(month), int(day), int(hour), int(min))        
        else: 
            raise ValueError()
        return dt        
    except ValueError:
        return None

def set_table(section5):
    max_level = struct.unpack_from('>H', section5, 15)[0]
    table = (
        -10, # define representative of level 0 (Missing Value)
        *struct.unpack_from('>'+str(max_level)+'H', section5, 18)
    )
    return np.array(table, dtype=np.int16)

def decode_runlength(code, hi_level):
    for raw in code:
        if raw <= hi_level:
            level = raw
            pwr = 0
            yield level
        else:
            length = (0xFF - hi_level)**pwr * (raw - (hi_level + 1))
            pwr += 1
            yield from repeat(level, length)

def load_jmara_grib2(file):
    with open(file, 'rb') as f:
        binary = f.read()

    len_ = {'sec0':16, 'sec1':21, 'sec3':72, 'sec4':82, 'sec6':6}

    end4 = len_['sec0'] + len_['sec1'] + len_['sec3'] + len_['sec4'] - 1
    len_['sec5'] = struct.unpack_from('>I', binary, end4+1)[0]
    section5 = binary[end4:(end4+len_['sec5']+1)]

    end6 = end4 + len_['sec5'] + len_['sec6']
    len_['sec7'] = struct.unpack_from('>I', binary, end6+1)[0]
    section7 = binary[end6:(end6+len_['sec7']+1)]

    highest_level = struct.unpack_from('>H', section5, 13)[0]
    level_table = set_table(section5)
    decoded = np.fromiter(decode_runlength(section7[6:], highest_level), dtype=np.int16).reshape((3360, 2560))

    transposed_flipped_data = np.flip(np.transpose(decoded), axis=1)
    return level_table[transposed_flipped_data] / 100
    
def read_msm(time):
    # 初期値から5時間後に配信される
    modeltime = time - offsets.Hour(5)
    # MSMは03シリーズ
    base_time = modeltime.replace(hour=modeltime.hour - (modeltime.hour % 3), minute=0, second=0)  
    # 対象時刻と初期値の時間差
    ft = (time - base_time).total_seconds() // 3600
    # 生存圏研究所ダウンロード元サイト
    http  = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original"  
    # データファイル名の指定
    day_dir = base_time.strftime("%Y/%m/%d")
    basename = "Z__C_RJTD_{}00_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin".format(base_time.strftime("%Y%m%d%H%M"))
    fname    = "./{}".format(basename)
    # すでにファイルが存在しなければ、ダウンロードを行う
    if os.path.exists(fname):
        print("{} Already exists.".format(basename))
    else:
        url      = "{}/{}/{}".format(http,  day_dir, basename)
        # wgetコマンドでデータのダウンロード
        subprocess.run("wget {} -P ./ > /dev/null 2>&1".format(url), shell=True)       
    grbs = pygrib.open(basename)
    # 海面気圧のデータを取得する
    prmsl_fc0 = grbs.select(parameterName='Pressure reduced to MSL', forecastTime=ft)[0]
    prmsl, lats, lons = prmsl_fc0.data()
    prmsl_flipped = np.flip(prmsl / 100, axis=0)
    # 気温のデータを取得する
    tmp_fc0 = grbs.select(parameterName='Temperature', forecastTime=ft)[0]
    tmp, lats, lons = tmp_fc0.data()
    tmp_flipped = np.flip(tmp - 273.15, axis=0)
    # 風速のデータを取得する
    u_fc0 = grbs.select(parameterName='u-component of wind', forecastTime=ft)[0]
    u, lats, lons = u_fc0.data()
    u_flipped = np.flip(u, axis=0)
    v_fc0 = grbs.select(parameterName='v-component of wind', forecastTime=ft)[0]
    v, lats, lons = v_fc0.data()
    v_flipped = np.flip(v, axis=0)
    return prmsl_flipped, tmp_flipped, u_flipped, v_flipped
    
# 緯度経度で指定したポイントの図上の座標などを取得する関数 
# 図法の座標 => pixel座標 => 図の座標　と3回の変換を行う
# pixel座標: plt.figureで指定した大きさxDPIに合わせ、左下を原点とするpixelで測った座標   
# 図の座標: axesで指定した範囲を(0,1)x(0,1)とする座標
# 3つの座標（図の座標, Pixel座標, 図法の座標）を出力する 

def transform_lonlat_to_figure(lonlat, ax, proj):
    # lonlat:経度と緯度  (lon, lat) 
    # ax: Axes図の座標系    例：fig.add_subplot()の戻り値
    # proj: axで指定した図法 
    #
    # 例：緯度経度をpointで与え、ステレオ図法にする場合
    # point = (140.0,35.0)
    # proj= ccrs.Stereographic(central_latitude=60, central_longitude=140) 
    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_subplot(1, 1, 1, projection=proj)
    # ax.set_extent([108, 156, 17, 55], ccrs.PlateCarree())
    #
    # 図法の変換：参照  https://scitools.org.uk/cartopy/docs/v0.14/crs/index.html                    
    point_proj = proj.transform_point(*lonlat, ccrs.PlateCarree())
    # pixel座標へ変換：参照　https://matplotlib.org/stable/tutorials/advanced/transforms_tutorial.html
    point_pix = ax.transData.transform(point_proj)
    # 図の座標へ変換                                                           
    point_fig = ax.transAxes.inverted().transform(point_pix)
    return point_fig, point_pix, point_proj

# 1地点のアメダスjsonデータから、elem要素で指定した値を返す(ただしFlagが 0以外は Noneとする)
# 要素:elem = 'temp','humidity','snow1h','snow6h','snow12h','snow24h','sun10m','sun1h','precipitation10m','precipitation1h','precipitation3h','precipitation24h','wind','windDirection'
def get_obs_value(amd_obs,elem):
    try:
        et = amd_obs[elem]
        if int(et[1]) != 0:
            return None
        return float(et[0])
    except Exception:
        return None
    
# 気象庁全国合成レーダーGPVのダウンロード
def download_time(time):    
    # 生存圏研究所ダウンロード元サイト
    http  = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/jma-radar/synthetic/original"  
    # データファイル名の指定
    day_dir = time.strftime("%Y/%m/%d")
    basename = "Z__C_RJTD_{}00_RDR_JMAGPV_Ggis1km_Prr05lv_ANAL_grib2.bin".format(time.strftime("%Y%m%d%H%M"))
    fname    = "./{}".format(basename)
    # すでにファイルが存在しなければ、ダウンロードを行う
    if os.path.exists(fname):
        print("{} Already exists.".format(basename))
    else:
        url      = "{}/{}/{}".format(http,  day_dir, basename)
        # wgetコマンドでデータのダウンロード
        subprocess.run("wget {} -P ./ > /dev/null 2>&1".format(url), shell=True)
        # ダウンロードした圧縮ファイルの解凍
        #subprocess.run("gunzip -xvf {} -C ./".format(fname), shell=True)   
    #GgisFile = "./Z__C_RJTD_{}00_RDR_GPV_Ggis0p25km_Pri60lv_Aper5min_ANAL_grib2.bin".format(time.strftime("%Y%m%d%H%M"))
    #return GgisFile
    return fname

# 描画指定：順に気圧(右上),気温(左上),湿球温度(右下),露点温度(左下))
npre_dispflag = False
temp_dispflag = False
wbt_dispflag = False
dp_dispflag = False

markersize_0 = 8 # マーカーサイズ
char_size=8 # 文字サイズ
barb_length=6 # 矢羽の長さ
dlon,dlat=1,1   # 緯線・経線は1度ごと

# 描画地域と描画時刻の設定
if len(sys.argv) == 2:
    dt = parse_datetime(arg)
else:
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)

# 描画開始メッセージ    
if dt:
    year=dt.year
    month=dt.month
    day=dt.day
    hour=dt.hour
    min=0
    dt = datetime(int(year), int(month), int(day), int(hour), int(min))
    print("読み込み観測時刻 {:4d}/{:02d}/{:02d} {:02d}:{:02d}".format(year,month,day,hour,min))
else:
    print('Usage: python script.py [YYYYMMDDHH(MM)]')
    exit()

# 対象のcouパラメータのリスト
cou_list = [2010, 2140, 2150, 2176, 2180, 2190, 2250, 5030, 5420]
 
# ヘッダー行を定義
header = ['ID', 'City', 'Country', 'Latitude', 'Longitude', 'Height', 'Temperature', 'Humidity', 'Pressure', 'Wind Direction', 'Wind Speed']
 
# CSVファイルにデータを書き込む
with open('weather_data.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)  # ヘッダー行を書き込む
   
    # 各URLについて処理を実行し、CSVに書き込む
    for cou in cou_list:
        data = process_url(cou, year, month, day, hour)
        for row in data:
            writer.writerow(row)

def read_hima(time, band):
  # Himawari-9
  # 日付とファイル名の生成
  time = time - offsets.Hour(9)
  day_dir = time.strftime("%Y%m/%d")
  basename = "NC_H09_{}_R21_FLDK.02401_02401.nc".format(time.strftime("%Y%m%d_%H%M"))

  # lftpコマンドを実行してFTPサーバーに接続
  PTree_ID = os.environ.get('PTree_ID')
  PTree_Pass = os.environ.get('PTree_Pass')

  # ダウンロードするファイルのURLを作成
  url = "ftp://ftp.ptree.jaxa.jp/jma/netcdf/{}/{}".format(day_dir, basename)

  # wgetコマンドを使用してファイルをダウンロード
  wget_command = "wget --user={} --password={} {} -P ./ > /dev/null 2>&1".format(PTree_ID, PTree_Pass, url)
  subprocess.run(wget_command, shell=True)

  # NetCDF ファイルを開く
  nc_file = nc.Dataset(basename, 'r')

  # 緯度、経度、およびデータの取得
  latitude = nc_file.variables['latitude'][:]
  longitude = nc_file.variables['longitude'][:]
  data = nc_file.variables[f"tbb_{band}"][:].reshape(2401, 2401)

  # メッシュグリッドを作成
  lon, lat = np.meshgrid(longitude, latitude)

  return data, lon, lat
    
# 観測データJSONの url作成
url_data_json= 'https://www.jma.go.jp/bosai/amedas/data/map/{:4d}{:02d}{:02d}{:02d}{:02d}00.json'
url_data_json=url_data_json.format(year,month,day,hour,min)
# 気象庁HPからアメダスデータを読み込む
response = urlopen(url_data_json)
content = response.read()
response.close()
data_json=content.decode()
dat_json = json.loads(data_json)

# アメダス地点Tableのurl
url_station_json="https://www.jma.go.jp/bosai/amedas/const/amedastable.json"
# アメダス地点Tableを読み込む
response = urlopen(url_station_json)
content = response.read()
response.close()
station_json=content.decode()
amd_json = json.loads(station_json)

# アメダスデータと同じ時刻のUTCを計算
time = pd.Timestamp(year,month,day,hour,min)
utc = time - offsets.Hour(9)

# 前1時間の雷実況
for i in range(1,12):
    time_liden = utc - offsets.Minute(5*i)

    # LIDENデータのURL
    data_url = "https://www.jma.go.jp/bosai/jmatile/data/nowc/{}00/none/{}00/surf/liden/data.geojson?id=liden"
    data_url=data_url.format(time_liden.strftime("%Y%m%d%H%M"),time_liden.strftime("%Y%m%d%H%M"))

    # データの取得
    response = requests.get(data_url)
    data = response.json()

    # データの解析
    lons_liden = []
    lats_liden = []

    for feature in data['features']:
        coordinates = feature['geometry']['coordinates']
        lon, lat = coordinates
        lons_liden.append(lon)
        lats_liden.append(lat)

### 解析雨量
# データのURL
data_url = "https://www.jma.go.jp/bosai/jmatile/data/rasrf/{}00/immed/{}00/surf/rasrf_point/data.geojson?id=rasrf_point"
data_url=data_url.format(utc.strftime("%Y%m%d%H%M"),utc.strftime("%Y%m%d%H%M"))

# データの取得
response = requests.get(data_url)
data = json.loads(response.text)
 
# 座標と値のリストを作成
coordinates = []
values = []
for feature in data["features"]:
    coordinate = feature["geometry"]["coordinates"]
    value = float(feature["properties"]["value"])
    coordinates.append(coordinate)
    values.append(value)
 
# 座標データをNumPy配列に変換
coordinates = np.array(coordinates)
x = coordinates[:, 0]
y = coordinates[:, 1]

# 値データをNumPy配列に変換
values = np.array(values)
 
# グリッドの作成
xi = np.linspace(np.min(x), np.max(x), 370)
yi = np.linspace(np.min(y), np.max(y), 481)
xi, yi = np.meshgrid(xi, yi)
 
# 値データを補間
zi = np.zeros_like(xi)
for i in range(len(values)):
    xi_index = np.abs(xi[0] - x[i]).argmin()
    yi_index = np.abs(yi[:, 0] - y[i]).argmin()
    zi[yi_index, xi_index] = values[i]
   
# GPVデータの時間の指定(年,月,日,時,分)
filepath = download_time(utc)

# データを読む
rain = load_jmara_grib2(filepath) #レーダーGPV
prmsl, tmp, u, v = read_msm(utc) #MSM海面気圧

# 地形データ取得
data = np.fromfile("LANDSEA.MSM_5K", dtype=np.dtype('>f4'))  # 地形バイナリデータをBig Endianの単精度浮動小数点数として読み込む

# データサイズに応じてreshape（MSMのグリッドサイズに合わせる）
data = data.reshape(505, 481)

# メッシュグリッドの作成
grid_lon_s, grid_lat_s = np.meshgrid(np.arange(120, 150 + 0.0625, 0.0625), np.arange(22.4, 47.6, 0.05))
sealand = np.flip(data*10000, axis=0)

#sealand[(grid_lon_s - grid_lat_s < 94.5) & (grid_lon_s < 132)] = 0
#sealand[(grid_lat_s > 45.5)] = 0
#sealand[(grid_lon_s > 145.5)] = 0

sealand = maximum_filter(sealand, size=(15, 15))

# ガウシアンフィルタを適用
sealand_filterd = gaussian_filter(sealand, sigma=1.0) # sigmaはガウス分布の標準偏差

# 図法指定                                                                             
proj = ccrs.PlateCarree()

# カラーバーの設定
#気象庁RGBカラー
jmacolors=np.array(
   [[0.95,0.95,0.95,1],#white
    [0.63,0.82,0.99,1],
    [0.13,0.55,0.99,1],
    [0.00,0.25,0.99,1],
    [0.98,0.96,0.00,1],
    [0.99,0.60,0.00,1],
    [0.99,0.16,0.00,1],
    [0.71,0.00,0.41,1]]
)
#等高線値
clevs = np.array([0,1,5,10,20,30,50,80]) 
   
## 入力ファイルの次元の指定（水平解像度約1km）
## 経度方向： 118°〜150° (1/80°間隔)
## 緯度方向： 20°〜48° (1/120°間隔)
mlon, mlat = 2560, 3360
slon, elon, rlon = 118, 150, 1/80
slat, elat, rlat = 20, 48, 1/120

for area in [0, 1, 2, 3, 4]:
    # 地図の描画範囲指定
    if (area == 0):
        i_area = [139, 147, 40, 46]
        areaname="Hokkaido"
    elif (area == 1):
        i_area = [134, 142, 33, 39]
        areaname="East"
    elif (area == 2):
        i_area = [128, 136, 31, 37]
        areaname="West"
    elif (area == 3):
        i_area = [135, 143, 36, 42]
        areaname="Tohoku"
    elif (area == 4):
        i_area = [123, 131, 25, 31]
        areaname="AmamiOkinawa"
        
    # 図のSIZE指定inch                                                                        
    fig = plt.figure(figsize=(8,6))
    # 余白設定                                                                                
    plt.subplots_adjust(left=0.04, right=1.1, bottom=0.0, top=1.0)                  
    # 作図                                                                                    
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.set_extent(i_area, proj)

    # カラーマップの作成
    norm = BoundaryNorm(clevs, len(clevs) - 1)

    # レーダーGPV描画
    lon = np.arange(slon, elon, rlon)
    lat = np.arange(slat, elat, rlat)
    LON, LAT = np.meshgrid(lon, lat)
    LON, LAT = LON.T, LAT.T
    #cs = ax.contourf(LON, LAT, rain, colors=jmacolors, levels=clevs, extend="max")
    cs = ax.contourf(xi, yi, zi, levels=clevs, cmap=ListedColormap(jmacolors), norm=norm)

    cb = plt.colorbar(cs, orientation="vertical", ticks=clevs, shrink=0.6)    
    #cb.ax.tick_params(labelsize=8) 

    # グリッド線を引く                                                               
    xticks=np.arange(-180,180,dlon)
    yticks=np.arange(-90,90.1,dlat)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)

    # 配列の宣言
    lat_list_t = []
    lon_list_t = []
    temp_list = []
    lat_list_p = []
    lon_list_p = []
    npre_list = []

    synop_file_name = "weather_data.csv"
    try:
        with open(synop_file_name, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                pressure = row['Pressure']
                if pressure != '-':
                    lat_list_p.append(float(row['Latitude']))
                    lon_list_p.append(float(row['Longitude']))
                    npre_list.append(float(pressure))
    except FileNotFoundError:
        print(f"File {synop_file_name} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # 地点プロット                                                                                                 
    for stno,val in dat_json.items():
        # 緯度・経度のtuple(度分形式)をtuple(度単位)に変換
        wlat = amd_json[stno]['lat'][0] + amd_json[stno]['lat'][1]/60.0
        wlon = amd_json[stno]['lon'][0] + amd_json[stno]['lon'][1]/60.0
        walt = amd_json[stno]['alt']
        # 天気
        weather = get_obs_value(val,'weather')
        # 風
        wind_ok = True
        ws = get_obs_value(val,'wind')
        wd = get_obs_value(val,'windDirection')
        if ws is not None and wd is not None:
            # 16方位風向と風速から、u,v作成   
            if ws == None or wd == None:
                u = None
                v = None
            else:
                wd = wd / 8.0 * math.pi  # 1/8 = (360/16) / 360 * 2
                au = -1.0 * ws * math.sin(wd)
                av = -1.0 * ws * math.cos(wd)
        else:
            wind_ok = False
        # 気温
        temp = get_obs_value(val,'temp')
        if temp is None:
            temp = np.nan
        elif walt < 800: #標高の高い観測点は無視する
            # 配列に格納
            tempsl = temp
            lat_list_t.append(wlat)
            lon_list_t.append(wlon)
            temp_list.append(tempsl)
        # 湿度
        hu = get_obs_value(val,'humidity')
        if hu is None:
            hu = -1.0
        # 露点温度
        dp_temp = -200.0
        if hu >= 0.0 and temp > -200.0:
            dp_temp = mpcalc.dewpoint_from_relative_humidity(temp * units.degC,hu/100.0).m
        # 更正気圧
        npre = get_obs_value(val,'normalPressure')
        if npre is None:
            npre = np.nan
        else:
            # 配列に格納
            lat_list_p.append(wlat)
            lon_list_p.append(wlon)
            npre_list.append(npre)   
        # 気圧
        pre = get_obs_value(val,'pressure')
        if pre is None:
            pre = -1.0

        # 湿球温度
        wb_temp = -200.0
        if dp_temp > -200.0 and temp > -200.0 and pre > 0.0:
            wb_temp = mpcalc.wet_bulb_temperature(pre * units.hPa, temp * units.degC, dp_temp * units.degC).m
        
        ## プロット
        fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj) 
        if ( fig_z[0] > 0.01 and fig_z[0] < 0.99  and fig_z[1] > 0.01 and fig_z[1] < 0.99):
            if weather == 0:
                color="orange"
            elif weather == 1:
                color="gray"
            elif weather == 7:
                color="green"
            elif weather == 10:
                color="blue"
            else:
                color="none"
            ax.plot(wlon, wlat, marker='o', markersize=markersize_0, color=color, transform=proj)
            if wind_ok and au*au+av*av>4.0: # 矢羽プロット
                ax.barbs(wlon, wlat, (au * units('m/s')).to('kt').m, (av * units('m/s')).to('kt').m, length=barb_length, transform=proj)
            if npre_dispflag and pre >= 0.0: # 気圧プロット
                ax.text(fig_z[0]+0.029, fig_z[1]+0.015,'{:6.1f}'.format(npre),size=char_size, color="black", transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
            if temp_dispflag and temp > -200.0: # 気温プロット
                ax.text(fig_z[0]-0.025, fig_z[1]+0.015,'{:5.1f}'.format(temp),size=char_size, color="red", transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
            if wbt_dispflag and wb_temp > -200.0: # 湿球温度プロット
                if dp_temp < 0:
                    color_temp = "purple"
                else:
                    color_temp = "black"
                ax.text(fig_z[0]+0.025, fig_z[1]-0.003,'{:5.1f}'.format(wb_temp),size=char_size, color=color_temp, transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")
            if dp_dispflag and dp_temp > -200.0: # 露点温度プロット
                if dp_temp < 0:
                    color_temp = "purple"
                else:
                    color_temp = "black"
                ax.text(fig_z[0]-0.025, fig_z[1]-0.003,'{:5.1f}'.format(dp_temp),size=char_size, color=color_temp, transform=ax.transAxes,verticalalignment="top", horizontalalignment="center")  

    prmsl = gaussian_filter(prmsl, sigma=4.0)
    tmp = gaussian_filter(tmp, sigma=4.0)

    # 線形補間
    grid_temp = griddata((lon_list_t, lat_list_t), temp_list, (grid_lon_s, grid_lat_s), method='linear')
    grid_npre = griddata((lon_list_p, lat_list_p), npre_list, (grid_lon_s, grid_lat_s), method='linear')
    grid_temp = np.where(np.isnan(grid_temp), tmp, grid_temp)
    grid_npre = np.where(np.isnan(grid_npre), prmsl, grid_npre)
    grid_npre = gaussian_filter(grid_npre, sigma=2.0)

    diff_npre = grid_npre - prmsl
    
    diff_npre = gaussian_filter(diff_npre, sigma=2.0)
    diff_npre[sealand_filterd > 1000.0] = grid_npre[sealand_filterd > 1000.0] - prmsl[sealand_filterd > 1000.0]
    diff_npre = gaussian_filter(diff_npre, sigma=2.0)

    grid_npre = prmsl + diff_npre
    
    #陸地から十分離れた格子は描画しない(MSMと実況の差が大きい場合があるため)
    #grid_npre[sealand_filterd <= 1] = np.nan
    #grid_temp[sealand_filterd <= 1] = np.nan

    # 描画領域のデータを切り出す（等圧線のラベルを表示するためのおまじない）
    lon_range = np.where((grid_lon_s[0, :] >= i_area[0] - 0.25) & (grid_lon_s[0, :] <= i_area[1] + 0.25))
    lat_range = np.where((grid_lat_s[:, 0] >= i_area[2] - 0.25) & (grid_lat_s[:, 0] <= i_area[3] + 0.25))

    # 切り出したい範囲のインデックスを取得
    lon_indices = lon_range[0]
    lat_indices = lat_range[0]

    # 切り出し
    grid_lon_sliced = grid_lon_s[lat_indices][:, lon_indices]
    grid_lat_sliced = grid_lat_s[lat_indices][:, lon_indices]
    psea_grid = grid_npre[lat_indices][:, lon_indices]
    temp_grid = grid_temp[lat_indices][:, lon_indices]

    # 等温線をプロット
    levels = np.arange(-30, 45, 3)
    #cont = plt.contour(grid_lon_sliced, grid_lat_sliced, temp_grid, levels=levels, linewidths=2, linestyles='solid', colors='red')

    # 等温線のラベルを付ける
    #plt.clabel(cont, fontsize=15)

    # 等圧線をプロット
    levels = np.arange(900, 1050, 1)
    cont = plt.contour(grid_lon_sliced, grid_lat_sliced, psea_grid, levels=levels, linewidths=2, colors='black')

    # 等圧線のラベルを付ける
    plt.clabel(cont, fontsize=15)

    # ベクトルの間引き間隔
    stride = 8
    ax.barbs(grid_lon_s[::stride, ::stride], grid_lat_s[::stride, ::stride], u[::stride, ::stride], v[::stride, ::stride], length=6, transform=proj)

    # LIDENプロット
    plt.scatter(lons_liden, lats_liden, marker='x', color='deeppink', s=200)

    ## H stamp
    maxid = detect_peaks(psea_grid, filter_size=40, dist_cut=10)
    for i in range(len(maxid[0])):
        wlon = grid_lon_sliced[0][maxid[1][i]]
        wlat = grid_lat_sliced[maxid[0][i]][0]
        # 図の範囲内に座標があるか確認                                                                           
        fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
        if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
            ax.plot(wlon, wlat, marker='x' , markersize=16, color="blue", transform=proj)
            ax.text(wlon - 0.12, wlat + 0.12, 'H', size=30, color="blue", transform=proj)
            val = psea_grid[maxid[0][i]][maxid[1][i]]
            ival = int(val)
            ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=24, color="blue", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

    ## L stamp
    minid = detect_peaks(psea_grid, filter_size=40, dist_cut=10, flag=1)
    for i in range(len(minid[0])):
        wlon = grid_lon_sliced[0][minid[1][i]]
        wlat = grid_lat_sliced[minid[0][i]][0]
        # 図の範囲内に座標があるか確認                                                                           
        fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
        if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
            ax.plot(wlon, wlat, marker='x' , markersize=16, color="red", transform=proj)
            ax.text(wlon - 0.12, wlat + 0.12, 'L', size=30, color="red", transform=proj)
            val = psea_grid[minid[0][i]][minid[1][i]]
            ival = int(val)
            ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=24, color="red", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

    # 海岸線
    ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
    # 図の説明
    plt.title('{}'.format("AMeDAS, RA1h, LIDEN1h"), loc='left',size=15)
    plt.title('Valid Time: {}'.format(dt), loc='right',size=15);
    #plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
    plt.savefig("latest{}.png".format(areaname), format="png")
    plt.clf()

char_size=16 # 文字サイズ
barb_length=8 # 矢羽の長さ
dlon,dlat=5,5   # 緯線・経線は1度ごと

# 描画地域と描画時刻の設定
if len(sys.argv) == 2:
    year = arg[:4]
    month = arg[4:6]
    day = arg[6:8]
    hour = arg[8:10]
else:
    jst = pytz.timezone('Asia/Tokyo')
    dt = datetime.now(jst) - timedelta(minutes=30)
    year=dt.year
    month=dt.month
    day=dt.day
    hour=dt.hour

dt = datetime(int(year), int(month), int(day), int(hour), 0)

# 初期値から5時間後に配信される
init = pd.Timestamp(year,month,day,hour,0) - offsets.Hour(6) 
init = init - offsets.Hour(init.hour % 3) 
ft = init + offsets.Hour(6)

# MSMは03シリーズ
base_time =  init - offsets.Hour(9) 

# データファイル名の指定
day_dir = base_time.strftime("%Y/%m/%d")
basename = "Z__C_RJTD_{}00_MSM_GPV_Rjp_L-pall_FH00-15_grib2.bin".format(base_time.strftime("%Y%m%d%H%M"))
url      = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/{}/{}".format(day_dir, basename)
subprocess.run("wget {} -P ./ > /dev/null 2>&1".format(url), shell=True)       
grbs = pygrib.open(basename)

# 配列の宣言
grid_lon_p, grid_lat_p = np.meshgrid(np.arange(120, 150 + 0.0625, 0.125), np.arange(22.4, 47.6, 0.1))

# データを取得する
height300 = gaussian_filter(np.flip(grbs.select(parameterName='Geopotential height', level=300, forecastTime=6)[0].data()[0], axis=0), sigma=4.0)
height500 = gaussian_filter(np.flip(grbs.select(parameterName='Geopotential height', level=500, forecastTime=6)[0].data()[0], axis=0), sigma=4.0)
tmp500 = np.flip(grbs.select(parameterName='Temperature', level=500, forecastTime=6)[0].data()[0] -273.15, axis=0)
u300 = np.flip(grbs.select(parameterName='u-component of wind', level=300, forecastTime=6)[0].data()[0], axis=0)
v300 = np.flip(grbs.select(parameterName='v-component of wind', level=300, forecastTime=6)[0].data()[0], axis=0)
u500 = np.flip(grbs.select(parameterName='u-component of wind', level=500, forecastTime=6)[0].data()[0], axis=0)
v500 = np.flip(grbs.select(parameterName='v-component of wind', level=500, forecastTime=6)[0].data()[0], axis=0)
u850 = np.flip(grbs.select(parameterName='u-component of wind', level=850, forecastTime=6)[0].data()[0], axis=0)
v850 = np.flip(grbs.select(parameterName='v-component of wind', level=850, forecastTime=6)[0].data()[0], axis=0)
u925 = gaussian_filter(np.flip(grbs.select(parameterName='u-component of wind', level=925, forecastTime=6)[0].data()[0], axis=0), sigma=4.0)
v925 = gaussian_filter(np.flip(grbs.select(parameterName='v-component of wind', level=925, forecastTime=6)[0].data()[0], axis=0), sigma=4.0)
# 相当温位
tmp700 = np.flip(grbs.select(parameterName='Temperature', level=700, forecastTime=6)[0].data()[0] -273.15, axis=0)
tmp850 = np.flip(grbs.select(parameterName='Temperature', level=850, forecastTime=6)[0].data()[0] -273.15, axis=0)
tmp925 = np.flip(grbs.select(parameterName='Temperature', level=925, forecastTime=6)[0].data()[0] -273.15, axis=0)
rh700 = np.flip(grbs.select(parameterName='Relative humidity', level=700, forecastTime=6)[0].data()[0], axis=0)
rh850 = np.flip(grbs.select(parameterName='Relative humidity', level=850, forecastTime=6)[0].data()[0], axis=0)
rh925 = np.flip(grbs.select(parameterName='Relative humidity', level=925, forecastTime=6)[0].data()[0], axis=0)
ept850 = gaussian_filter(mpcalc.equivalent_potential_temperature(850*units('hPa'), (tmp850+273.15) * units('K'), mpcalc.dewpoint_from_relative_humidity((tmp850+273.15) * units('K'), rh850 / 100)), sigma=4.0)
ept925 = gaussian_filter(mpcalc.equivalent_potential_temperature(925*units('hPa'), (tmp925+273.15) * units('K'), mpcalc.dewpoint_from_relative_humidity((tmp925+273.15) * units('K'), rh850 / 100)), sigma=4.0)
ttd = (tmp700 - mpcalc.dewpoint_from_relative_humidity((tmp700+273.15) * units('K'), rh700 / 100).magnitude)
# kindex = tmp850 - tmp500 + mpcalc.dewpoint_from_relative_humidity((tmp850+273.15) * units('K'), rh850 / 100).magnitude - ttd
kindex = 20 - tmp500 + mpcalc.dewpoint_from_relative_humidity((tmp850+273.15) * units('K'), rh850 / 100).magnitude - ttd

tmp500 = gaussian_filter(tmp500, sigma=4)
tmp850 = gaussian_filter(tmp850, sigma=4)

data_wv, lon, lat = read_hima(ft, '08')
data_ir, lon, lat = read_hima(time, '13')

# 図法指定                                                                             
proj = ccrs.PlateCarree()                                                      
# 図のSIZE指定inch                                                                        
fig = plt.figure(figsize=(8,6))
# 余白設定                                                                                 
plt.subplots_adjust(left=0.14, right=0.9, bottom=0.0, top=0.9) 

### 計算 ###
lats = np.arange(22.4, 47.6, 0.1)
lons = np.arange(120, 150 + 0.0625, 0.125)

dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
ug300, vg300 = mpcalc.geostrophic_wind(height300 * units('m'), dx=dx, dy=dy, latitude=grid_lat_p * units('degree'))
ug500, vg500 = mpcalc.geostrophic_wind(height500 * units('m'), dx=dx, dy=dy, latitude=grid_lat_p * units('degree'))

vor = mpcalc.vorticity(u500 * units('m/s'), v500 * units('m/s'), dx=dx, dy=dy) * 1000000
vor = gaussian_filter(vor, sigma=2.0)

### 300hPa ###
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
cont = plt.contour(grid_lon_p, grid_lat_p, height300, levels=np.arange(5100, 6000, 60), linewidths=2, colors='black')
plt.clabel(cont, fontsize=15)

# 描画
data = np.flipud(data_wv)
plt.imshow(data_wv, cmap='gray_r', extent=(lon.min(), lon.max(), lat.max(), lat.min()), origin='lower', transform=proj)
plt.tight_layout(rect=[0, 0, 1, 0.96])

# 風速の計算
#wind_speed = mpcalc.wind_speed(ug300 * units('m/s'), vg300 * units('m/s')).to(units.knots)
wind_speed = mpcalc.wind_speed(u300 * units('m/s'), v300 * units('m/s')).to(units.knots)
wind_speed = gaussian_filter(wind_speed, sigma=4.0)
#plt.contourf(grid_lon_p, grid_lat_p, wind_speed, levels=[0, 60, 120, np.inf], colors=['none', 'blue', 'purple'], alpha=0.2)
cont = plt.contour(grid_lon_p, grid_lat_p, wind_speed, levels=[60], linewidths=2, colors='lightblue')
cont = plt.contour(grid_lon_p, grid_lat_p, wind_speed, levels=[80], linewidths=2, colors='lightpink')
cont = plt.contour(grid_lon_p, grid_lat_p, wind_speed, levels=[100], linewidths=2, colors='violet')

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("WindSpeed300, WV Image"), loc='left',size=15)
plt.title('Valid Time: {}'.format(ft), loc='right',size=15);
#plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
plt.savefig("latest_300.png", format="png")
plt.clf()

### 500hPa ###
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
cont = plt.contour(grid_lon_p, grid_lat_p, tmp500, levels=np.arange(-60, 60, 3), linewidths=2, linestyles='solid', colors='red')
plt.clabel(cont, fontsize=15)
cont = plt.contour(grid_lon_p, grid_lat_p, height500, levels=np.arange(5100, 6000, 60), linewidths=2, colors='black')
plt.clabel(cont, fontsize=15)

plt.contourf(grid_lon_p, grid_lat_p, vor, levels=[-float('inf'), 0, 100, float('inf')], colors=['none', 'darkorange', 'brown'])

# ベクトルの間引き間隔
stride = 10
ax.barbs(grid_lon_p[::stride, ::stride], grid_lat_p[::stride, ::stride], u500[::stride, ::stride], v500[::stride, ::stride], length=4, transform=proj)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("Z500, T500, Vor500"), loc='left',size=15)
plt.title('Valid Time: {}'.format(ft), loc='right',size=15);
#plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
plt.savefig("latest_500.png", format="png")
plt.clf()

### 700hPa ###
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
cont = plt.contour(grid_lon_p, grid_lat_p, tmp850, levels=np.arange(-60, 60, 3), linewidths=2, linestyles='solid', colors='red')
plt.clabel(cont, fontsize=15)

plt.contourf(grid_lon_p, grid_lat_p, ttd, levels=[-float('inf'), 3, 15, float('inf')], colors=['lightgreen', 'none', 'yellow'])

# ベクトルの間引き間隔
stride = 10
ax.barbs(grid_lon_p[::stride, ::stride], grid_lat_p[::stride, ::stride], u850[::stride, ::stride], v850[::stride, ::stride], length=4, transform=proj)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("TTd700, T850, Wind850"), loc='left',size=15)
plt.title('Valid Time: {}'.format(ft), loc='right',size=15);
#plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
plt.savefig("latest_700.png", format="png")
plt.clf()

### 850hPa ###
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
cont = plt.contour(grid_lon_p, grid_lat_p, ept850, levels=np.arange(210, 390, 3), linewidths=1, linestyles='solid', colors='green')
cont2 = plt.contour(grid_lon_p, grid_lat_p, ept850, levels=np.arange(210, 390, 15), linewidths=2, linestyles='solid', colors='green')
plt.clabel(cont, fontsize=15)

plt.contourf(grid_lon_p, grid_lat_p, kindex, levels=[-float('inf'), 10, 25, 40, np.inf], colors=['none', 'yellow', 'pink', 'red'])

# ベクトルの間引き間隔
stride = 10

ax.barbs(grid_lon_p[::stride, ::stride], grid_lat_p[::stride, ::stride], u850[::stride, ::stride], v850[::stride, ::stride], length=4, transform=proj)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("EPT850, Wind850, K-index"), loc='left',size=15)
plt.title('Valid Time: {}'.format(ft), loc='right',size=15);
#plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
plt.savefig("latest_850.png", format="png")
plt.clf()

### Surface ###
# 作図                                                             
plt.subplots_adjust(left=0.04, right=1.1, bottom=0.0, top=0.96)  
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
prmsl = gaussian_filter(prmsl, sigma=8.0) 
cont = plt.contour(grid_lon_s, grid_lat_s, prmsl, levels=np.arange(900, 1100, 4), linewidths=2, linestyles='solid', colors='pink')
plt.clabel(cont, fontsize=15)

data = np.flipud(data_ir)
plt.imshow(data_ir, cmap='gray_r', extent=(lon.min(), lon.max(), lat.max(), lat.min()), origin='lower', transform=proj)

cs = ax.contourf(LON, LAT, rain, colors=jmacolors, levels=clevs, extend="max")
cb = plt.colorbar(cs, orientation="vertical", ticks=clevs, shrink=0.6)    
cb.ax.tick_params(labelsize=8)

## H stamp
maxid = detect_peaks(prmsl, filter_size=80, dist_cut=30)
for i in range(len(maxid[0])):
    wlon = grid_lon_s[0][maxid[1][i]]
    wlat = grid_lat_s[maxid[0][i]][0]
    # 図の範囲内に座標があるか確認                                                                           
    fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
    if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
        ax.plot(wlon, wlat, marker='x' , markersize=16, color="blue", transform=proj)
        ax.text(wlon - 0.6, wlat + 0.6, 'H', size=30, color="blue", transform=proj)
        val = prmsl[maxid[0][i]][maxid[1][i]]
        ival = int(val)
        ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=24, color="blue", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

## L stamp
minid = detect_peaks(prmsl, filter_size=80, dist_cut=30, flag=1)
for i in range(len(minid[0])):
    wlon = grid_lon_s[0][minid[1][i]]
    wlat = grid_lat_s[minid[0][i]][0]
    # 図の範囲内に座標があるか確認                                                                           
    fig_z, _, _ = transform_lonlat_to_figure((wlon,wlat),ax,proj)
    if ( fig_z[0] > 0.05 and fig_z[0] < 0.95  and fig_z[1] > 0.05 and fig_z[1] < 0.95):
        ax.plot(wlon, wlat, marker='x' , markersize=16, color="red", transform=proj)
        ax.text(wlon - 0.6, wlat + 0.6, 'L', size=30, color="red", transform=proj)
        val = prmsl[minid[0][i]][minid[1][i]]
        ival = int(val)
        ax.text(fig_z[0], fig_z[1] - 0.025, str(ival), size=24, color="red", transform=ax.transAxes, verticalalignment="top", horizontalalignment="center")

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("Psea, Rader and IR Image"), loc='left',size=15)
plt.title('Valid Time: {}'.format(time), loc='right',size=15);
#plt.savefig("{}.png".format(time.strftime("%Y%m%d%H%M")), format="png")
plt.savefig("latest.png", format="png")
plt.clf()
