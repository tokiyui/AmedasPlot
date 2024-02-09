#!/usr/bin/env python
# coding: utf-8

import pygrib, pytz, subprocess, sys, os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import pandas.tseries.offsets as offsets
from datetime import datetime, timedelta
from metpy.units import units
from scipy.ndimage import gaussian_filter
import netCDF4 as nc
from ftplib import FTP

markersize_0 = 2 # マーカーサイズ
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
init = pd.Timestamp(year,month,day,hour,0) - offsets.Hour(5) 
init = init - offsets.Hour(init.hour % 3) 
ft = init + offsets.Hour(6)
# MSMは03シリーズ
base_time =  init - offsets.Hour(9) 

# データファイル名の指定
day_dir = base_time.strftime("%Y/%m/%d")
basename = "Z__C_RJTD_{}00_MSM_GPV_Rjp_L-pall_FH00-15_grib2.bin".format(base_time.strftime("%Y%m%d%H%M"))
url      = "http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/{}/{}".format(day_dir, basename)
subprocess.run("wget -q {} -P ./".format(url), shell=True)       
grbs = pygrib.open(basename)

# 配列の宣言
grid_lon, grid_lat = np.meshgrid(np.arange(120, 150 + 0.0625, 0.125), np.arange(22.4, 47.6, 0.1))

# データを取得する
height = np.flip(grbs.select(parameterName='Geopotential height', level=500, forecastTime=6)[0].data()[0], axis=0)
tmp = np.flip(grbs.select(parameterName='Temperature', level=500, forecastTime=6)[0].data()[0] - 273.15, axis=0)
u = np.flip(grbs.select(parameterName='u-component of wind', level=850, forecastTime=6)[0].data()[0], axis=0)
v = np.flip(grbs.select(parameterName='v-component of wind', level=850, forecastTime=6)[0].data()[0], axis=0)
# 相当温位
tmp850 = np.flip(grbs.select(parameterName='Temperature', level=850, forecastTime=6)[0].data()[0] * units('K'), axis=0)
rh = np.flip(grbs.select(parameterName='Relative humidity', level=850, forecastTime=6)[0].data()[0] / 100, axis=0)
dewpoint = mpcalc.dewpoint_from_relative_humidity(tmp850, rh)
ept = mpcalc.equivalent_potential_temperature(850*units('hPa'), tmp850, dewpoint)

height = gaussian_filter(height, sigma=4.0)
tmp = gaussian_filter(tmp, sigma=4.0)
ept = gaussian_filter(ept, sigma=4.0)

# Himawari-9
# 日付とファイル名の生成
day_dir = base_time.strftime("%Y%m/%d")
basename = "NC_H09_{}_R21_FLDK.02401_02401.nc".format(base_time.strftime("%Y%m%d_%H%M"))

# lftpコマンドを実行してFTPサーバーに接続
PTree_ID = os.environ.get('PTree_ID')
PTree_Pass = os.environ.get('PTree_Pass')

# ダウンロードするファイルのURLを作成
url = "ftp://ftp.ptree.jaxa.jp/jma/netcdf/{}/{}".format(day_dir, basename)

# wgetコマンドを使用してファイルをダウンロード
wget_command = "wget -q --user={} --password={} {} -P ./".format(PTree_ID, PTree_Pass, url)
subprocess.run(wget_command, shell=True)

# NetCDF ファイルを開く
nc_file = nc.Dataset(basename, 'r')

# 緯度、経度、およびデータの取得
latitude = nc_file.variables['latitude'][:]
longitude = nc_file.variables['longitude'][:]
data = nc_file.variables['tbb_13'][:].reshape(2401, 2401)

# メッシュグリッドを作成
lon, lat = np.meshgrid(longitude, latitude)

# ファイルを閉じる
nc_file.close()

# データをサンプリングする
sample = 5
sampled_data = data[::sample, ::sample]  # 2つおきにサンプリング

# サンプリングされたデータに対応する緯度経度を抽出する
sampled_lon = lon[::sample, ::sample]
sampled_lat = lat[::sample, ::sample]


# 図法指定                                                                             
proj = ccrs.PlateCarree()
# 図のSIZE指定inch                                                                        
fig = plt.figure(figsize=(8,6))
# 余白設定                                                                                
plt.subplots_adjust(left=0.04, right=1.05, bottom=0.0, top=1.0)                  
# 作図                                                                                    
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([120, 150, 22.4, 47.6], proj)

# グリッド線を引く                                                               
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, alpha=0.8)
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,5))
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,5))

# プロット
cont = plt.contour(grid_lon, grid_lat, ept, levels=np.arange(210, 390, 9), linewidths=2, linestyles='solid', colors='green')
plt.clabel(cont, fontsize=15)

cont = plt.contour(grid_lon, grid_lat, height, levels=np.arange(5100, 6000, 60), linewidths=2, colors='black')
plt.clabel(cont, fontsize=15)

#plt.contourf(grid_lon, grid_lat, tmp, cmap='turbo', levels=np.arange(-48, 3, 3))
#cb = plt.colorbar(orientation="vertical", shrink=0.6)    
#cb.ax.tick_params(labelsize=8)

print(data)
# 描画
#plt.contourf(lon, lat, data, cmap='gray_r')
plt.contourf(sampled_lon, sampled_lat, sampled_data, cmap='gray_r')
#plt.title('Brightness Temperature - Band tbb_08')

# ベクトルの間引き間隔
stride = 5

# データを間引く
grid_lon_sparse = grid_lon[::stride, ::stride]
grid_lat_sparse = grid_lat[::stride, ::stride]
u_sparse = u[::stride, ::stride]
v_sparse = v[::stride, ::stride]

ax.barbs(grid_lon_sparse, grid_lat_sparse, u_sparse, v_sparse, length=4, transform=proj)

# 海岸線
ax.coastlines(resolution='10m', linewidth=1.6, color='black')  
            
# 図の説明
plt.title('{}'.format("Z500, T500, EPT850, UV850"), loc='left',size=15)
plt.title('Valid Time: {}'.format(ft), loc='right',size=15);
#plt.savefig("{}.jpg".format(time.strftime("%Y%m%d%H%M")), format="jpg")
plt.savefig("latest_upper.jpg", format="jpg")
