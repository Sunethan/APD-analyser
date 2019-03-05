# APD-analyser

## 目的
這是專門用來模擬作為單光子偵測器之雪崩光電二極體元件的程式。因為目前碩班有在使用 Sentaurus TCAD 半導體元件模擬軟體，模擬 APD 的結構及其製程，所以為了確保該模擬符合我們在論文或課本上學到的 Avalanche Breakdown 現象，我開始用獨立於模擬軟體的方式——使用 Python 運算模擬結構的崩潰電壓——來確認模擬結果是否可靠，後來兩者結果是一致符合的，崩潰電壓大多誤差落於 0~5(V) 左右，其隨濃度變化趨勢都相同。而在算了許許多多次之後，決定寫為一個可套用至任何摻質濃度的 GUI，這樣就不需要每次為了新的濃度分佈來重寫程式碼——總是反覆地讀取CSV、建陣列存資料、計算結構參數、空乏區寬度與崩潰電壓。話說回來，這個程式目前預設的結構為 InP/InGaAsP/InGaAs 的 SAM-APD 的一維結構，裡頭的參數皆參考自諸多論文，目前尚未有時間整理這些模擬參數的由來，有興趣了解者可以再跟我聯繫。因為時間有限，所以還有功能尚未很齊全，但基本上可以做如下模擬：
這是專門用來模擬作為單光子偵測器之雪崩光電二極體元件的程式。因為目前碩班有在使用 Sentaurus TCAD 半導體元件模擬軟體，模擬 APD 的結構及其製程，所以為了確保該模擬符合我們在論文或課本上學到的 Avalanche Breakdown 現象，我開始用獨立於模擬軟體的方式——使用 Python 運算模擬結構的崩潰電壓，後來確定兩者結果是一致符合的。而在算了許許多多次之後，決定寫為一個可套用至任何摻質濃度的 GUI，這樣就不需要每次為了新的濃度分佈來重寫程式碼——總是反覆地讀取CSV、建陣列存資料、計算結構參數、空乏區寬度與崩潰電壓。話說回來，這個程式目前預設的結構為 InP/InGaAsP/InGaAs 的 SAM-APD 的一維結構，裡頭的參數皆參考自諸多論文，目前尚未有時間整理這些模擬參數的由來，有興趣了解者可以再跟我聯繫。因為時間有限，所以還有功能尚未很齊全，但基本上可以做如下模擬：

## 使用畫面
![image](https://github.com/Sunethan/APD-analyser/blob/master/readme-figures/Doping-Profile-fitting.png)
![image](https://github.com/Sunethan/APD-analyser/blob/master/readme-figures/depletion-width.png)
![image](https://github.com/Sunethan/APD-analyser/blob/master/readme-figures/Statistics.jpg)

## 功能
1. 輸入任意的摻質濃度分佈。
> 這是用在 InP / InGaAsP / InGaAs 的 SAM-APD 結構，所以雖然濃度分佈是任意的，但必須指名清楚的化合物分佈範圍。例如說，在 x < -3.6 與 x > -0.45 的範圍就是 InP。只要將範圍設定好，那麼模擬時用到的物理參數就是正確的。

2. 設定新網格。
> 因為摻質濃度分佈數據點可能不夠密集，或者是太密集，進而使得後續模擬有顯著偏差，所以可以藉由給定的網格線（mesh line），藉由線性近似來設定新的一組濃度分佈數據。

3. 藉由空乏近似（Depletion approximation），計算空乏區寬度隨偏壓的變化
> 繪製空乏區寬度隨偏壓的原因是用來優化在尋找崩潰電壓時所需的網格分佈。這部分比較複雜，未來如有機會可再詳細說明。倘若沒有優化好，那麼在尋找崩潰電壓時會很容易不收斂。

4. 尋找發生雪崩崩潰的電壓值（avalanche breakdown voltage）
> 基本上是透過倍增積分（multiplication integral）是否等於一來判斷是否已達崩潰。

## 結構
* main.py：主要的 tkinter 執行檔
* physfunc.py：所有物理模型與相關計算的檔案
* defaults.csv：開啟程式時所載入的預設值，主要是最近一次模擬時於各欄位（tk.entry）中的數值。
* VbLog.csv：「顯示崩潰電壓運算結果」的暫存檔
* VptLog.csv：「顯示擊穿電壓運算結果」的暫存檔
* log/Allsaved.csv：「每次模擬紀錄的各欄位數值」路徑檔
* log/Statistics.csv：顯示於統計列表的內容
* log/TempStatistics.csv：更新統計列表之暫存檔
* log/2019-mm-dd-hh-mm-ss-Settings.csv：該次模擬紀錄檔
* Sample_DopingProfile/DopingProfile_Total.csv：摻質濃度分佈範例檔

## 尚未實現之功能

* 還沒有做好使用者引導介面〈User Guide〉
* 因為還沒有設計藉由輸入之二維摻質濃度並自動求出曲率半徑的演算法，所以〈Breakdown voltage〉的 Cylindrical 按鈕沒有效果。
* 無法不透過重新計算擊穿或崩潰電壓來更新〈Statstics〉中的統計資料列表。
* 擊穿 InP/InGaAsP 界面與 InGaAsP/InGaAs 界面之擊穿電壓並不相同，雖然有顯示於左下之〈Variable History〉欄位中（Vpt 與 Vpt2），但並不會將 Vpt2 更新於統計列表中。
