import numpy as np
def binarize_by_ratio(arr,ratio):
    #acording to ratio binarize the array
    #import numpy as np
    # 创建一个示例二维数组
    #arr = np.random.rand(4, 5)  # 生成一个4x5的随机数组
    # 打印原始数组
    # print("Original array:")
    # print(arr)
    # 将二维数组展平成一维数组

    flat_arr = arr.flatten()
    # 对一维数组进行排序并获取排序后的索引
    sorted_indices = np.argsort(flat_arr)
    # 计算前30%和后70%的分界索引
    # front ratio to set 0
    cutoff = int(len(flat_arr) * ratio)
    # 创建一个新的一维数组
    new_flat_arr = np.ones_like(flat_arr)
    # 将前30%对应的索引位置置为0
    new_flat_arr[sorted_indices[:cutoff]] = 0
    # 将一维数组重新变回二维数组
    new_arr = new_flat_arr.reshape(arr.shape)
    # 打印处理后的数组
    #print("Processed array:")
    #print(new_arr)
    return(new_arr)