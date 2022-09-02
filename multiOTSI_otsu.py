from cmath import nan
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal
import cv2
from scipy.cluster.hierarchy import dendrogram , linkage
from skimage.filters import threshold_multiotsu
import argparse
import shapefile 
from sqlalchemy.sql.expression import true
import os
from matplotlib.pyplot import MultipleLocator, thetagrids
from sklearn.cluster import KMeans,AgglomerativeClustering
import math
import copy
def euler_distance(point1: np.ndarray, point2: list) -> float:
    """
    计算两点之间的欧拉距离，支持多维
    """
    distance = 0.0
    for a, b in zip(point1, point2):
        distance += math.pow(a - b, 2)
    return math.sqrt(distance)

class ClusterNode(object):
    def __init__(self, vec, left=None, right=None, distance=-1, id=None, count=1):
        """
        :param vec: 保存两个数据聚类后形成新的中心
        :param left: 左节点
        :param right:  右节点
        :param distance: 两个节点的距离
        :param id: 用来标记哪些节点是计算过的
        :param count: 这个节点的叶子节点个数
        """
        self.vec = vec
        self.left = left
        self.right = right
        self.distance = distance
        self.id = id
        self.count = count

class Hierarchical(object):
    def __init__(self, k = 1):
        assert k > 0
        self.k = k
        self.labels = None
    def fit(self, x):
        nodes = [ClusterNode(vec=v, id=i) for i,v in enumerate(x)]
        distances = {}
        point_num, future_num = np.shape(x)  # 特征的维度
        self.labels = [ -1 ] * point_num
        currentclustid = -1
        while len(nodes) > self.k:
            min_dist = math.inf
            nodes_len = len(nodes)
            closest_part = None  # 表示最相似的两个聚类
            for i in range(nodes_len - 1):
                for j in range(i + 1, nodes_len):
                    # 为了不重复计算距离，保存在字典内
                    d_key = (nodes[i].id, nodes[j].id)
                    if d_key not in distances:
                        distances[d_key] = euler_distance(nodes[i].vec, nodes[j].vec)
                    d = distances[d_key]
                    if d < min_dist:
                        min_dist = d
                        closest_part = (i, j)
            # 合并两个聚类
            part1, part2 = closest_part
            node1, node2 = nodes[part1], nodes[part2]
            new_vec = [ (node1.vec[i] * node1.count + node2.vec[i] * node2.count ) / (node1.count + node2.count)
                        for i in range(future_num)]
            new_node = ClusterNode(vec=new_vec,
                                   left=node1,
                                   right=node2,
                                   distance=min_dist,
                                   id=currentclustid,
                                   count=node1.count + node2.count)
            currentclustid -= 1
            del nodes[part2], nodes[part1]   # 一定要先del索引较大的
            nodes.append(new_node)
        self.nodes = nodes
        self.calc_label()

    def calc_label(self):
        """
        调取聚类的结果
        """
        for i, node in enumerate(self.nodes):
            # 将节点的所有叶子节点都分类
            self.leaf_traversal(node, i)

    def leaf_traversal(self, node: ClusterNode, label):
        """
        递归遍历叶子节点
        """
        if node.left == None and node.right == None:
            self.labels[node.id] = label
        if node.left:
            self.leaf_traversal(node.left, label)
        if node.right:
            self.leaf_traversal(node.right, label)

def read_img(filename):
    dataset = gdal.Open(filename)  # 打开文件
    im_width = dataset.RasterXSize  # 栅格矩阵的列数
    im_height = dataset.RasterYSize  # 栅格矩阵的行数
    im_geotrans = dataset.GetGeoTransform()  # 仿射矩阵
    im_proj = dataset.GetProjection()  # 地图投影信息
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
    #im_data[im_data<0] = 0
    #im_data , precentage = stretch_n(im_data , 0 , 255)
    #im_data = im_data.astype(np.uint8)  # 将数据写成数组，对应栅格矩阵
    #im_data[im_data ==0 ] = -1
    del dataset  # 关闭对象，文件dataset
    #return im_proj, im_geotrans, im_data, im_height, im_width, precentage
    return im_proj, im_geotrans, im_data, im_height, im_width

def stretch_n(bands, img_min, img_max, lower_percent=0, higher_percent=100):
    """
    :param bands:  目标数据，numpy格式
    :param img_min:   目标位深的最小值，以8bit为例，最大值为255， 最小值为0
    :param img_max:    目标位深的最大值
    :return:
    """
    out = np.zeros_like(bands).astype(np.float32)
    a = img_min
    b = img_max
    c = np.percentile(bands[:, :], lower_percent)
    d = np.percentile(bands[:, :], higher_percent)
    t = a + (bands[:, :] - c) * (b - a) / (d - c)
    t[t < a] = a
    t[t > b] = b
    out[:, :] = t
    precentage = (d-c) / (b-a)
    return out, precentage

def write_img(filename, im_proj, im_geotrans, im_data):

    # gdal数据类型包括
    # gdal.GDT_Byte,
    # gdal .GDT_UInt16, gdal.GDT_Int16, gdal.GDT_UInt32, gdal.GDT_Int32,
    # gdal.GDT_Float32, gdal.GDT_Float64

    # 判断栅格数据的数据类型
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    # 判读数组维数
    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    else:
        im_bands, (im_height, im_width) = 1, im_data.shape

    # 创建文件
    driver = gdal.GetDriverByName("GTiff")  # 数据类型必须有，因为要计算需要多大内存空间
    dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)
    if(dataset!= None):
        dataset.SetGeoTransform(im_geotrans) #写入仿射变换参数
        dataset.SetProjection(im_proj)  
    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)  # 写入数组数据
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
    del dataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', nargs='+', type=str, default=r"E:\14_HHK\YRDS2SSSINOwater.tif", help='img path(s)')
    parser.add_argument('--output', nargs='+', type=str, default=r"E:\14_HHK\YRDS2SSSINOwaterRESULT.tif", help='img save path(s)')
    parser.add_argument('--visualize', type=str, default= True, help='visualize_result')  # file/folder, 0 for webcam
    opt = parser.parse_args()
    # Setting the font size for all plots.
    matplotlib.rcParams['font.size'] = 9
    #path = r"D:\huanghekou\10_34Y_NDVIbuffer\1GEEdownload_NDVImin\NNMB_mNDWI\15NNMB.tif"
    # The input image.
    print('_____________开始读取数据_______________')
    proj, geotrans, image, row, column  = read_img(opt.source)#读取数据
    # Applying multi-Otsu threshold for the default value, generating
    # three classes.
    print('_____________开始分类数据&保存_______________')
    image_use = image[image>-1]
    thresholds = threshold_multiotsu(image_use,classes=3)
    print("大津算法的阈值为{}".format(thresholds))
    index_image = np.where((image>-1)==True)
    #strech_thresholds = thresholds*precentage
    #print('0-255拉伸阈值{}'.format(thresholds))
    #print('NNMB阈值{}'.format(strech_thresholds))
    # Using the threshold values, we generate the three regions.
    image_f = image_use.copy().flatten()
    regions = np.digitize(image_f, bins=thresholds)#设置了阈值  np.array([6.81])
    image_regions = np.full(image.shape,np.NAN)
    for num , (i , j) in enumerate(zip(index_image[0],index_image[1])):
        image_regions[i][j] = regions[num]
    write_img(opt.output ,proj, geotrans ,image_regions.astype(np.uint8))#保存数据
    #print(regions)#打印阈值
    print('_____________可视化数据_______________')
    if opt.visualize :
        fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(40, 3))

        # Plotting the original image.
        ax[0].imshow(image, cmap='gray')
        ax[0].set_title('Original')
        ax[0].axis('off')

        #x_major_locator=MultipleLocator(0.1)
        #把x轴的刻度间隔设置为1，并存在变量里
        #y_major_locator=MultipleLocator(10)

        # Plotting the histogram and the two thresholds obtained from
        # multi-Otsu.
        ax[1].hist((image.ravel()[image.ravel()>-999]),bins=100)
        ax[1].set_title('Histogram')
        #ax[1].xaxis.set_major_locator(x_major_locator)
        for thresh in thresholds:
            ax[1].axvline(thresh, color='r')
        image_clip_nodata = image.ravel()[-999<image.ravel()]
        ax[2].hist(image_clip_nodata.ravel()[image_clip_nodata.ravel()<10], bins=100)
        ax[2].set_title('Histogram_partview')
        #ax[2].xaxis.set_major_locator(x_major_locator)
        for thresh in thresholds:
            ax[2].axvline(thresh, color='r')
        # Plotting the Multi Otsu result.
        ax[3].imshow(image_regions, cmap='jet')
        ax[3].set_title('Multi-Otsu result')
        ax[3].axis('off')

        plt.subplots_adjust()

        plt.show()