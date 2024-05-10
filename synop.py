import requests
import re
import csv
import sys
 
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
        pressure = match[9].replace(' Hpa', '') if 'N/A' not in match[9] else 'N/A'
       
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
 
# 引数から年月日時を取得
if len(sys.argv) != 2:
    print("Usage: python script.py YYYYMMDDHH")
    sys.exit(1)

datetime_arg = sys.argv[1]
year = datetime_arg[:4]
month = datetime_arg[4:6]
day = datetime_arg[6:8]
hour = datetime_arg[8:]

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
