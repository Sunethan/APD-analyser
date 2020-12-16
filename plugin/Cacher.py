# 快取器類別
class Cacher:

    def __init__(self):
        self.computed = False
        self.data = None

    def save(self, data):
        # 把資料存到記憶體中，這邊你也可以自行實作把資料存到硬碟
        # !目前版本不會檢查輸入，理論上輸入不同應該要重新算一次
        print("save data to cache")
        self.data = data
        self.computed = True

    def read(self):
        # 從記憶體中讀取資料
        print("read data from cache")
        return self.data

    @property
    def is_computed(self):
        # 回傳記憶體目前是否有算好的資料
        return self.computed

    def clean_cache(self):
        # 清除記憶體快取資料
        self.data = None
        self.computed = False


# 快取裝飾器
def cached(cacher: Cacher):
    def decorator(func):
        def cached_func(*arg, **kwarg):

            if cacher.is_computed:
                # 如果有算過就從快取拿資料
                result = cacher.read()
                return result
            else:
                # 如果沒算過就從頭算，然後把算完的結果存快取
                result = func(*arg, **kwarg)
                cacher.save(result)
                return result

        return cached_func

    return decorator


# 建立快取器物件實體
cacher = Cacher()


# 用快取物件封裝add_one函數
@cached(cacher)
def add_one(x):
    print("compute y=x+1")
    y = x + 1
    return y