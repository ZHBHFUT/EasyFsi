# EasyFsi python interface

[toc]

## 1. Constants

### `FieldLocation`

变量存储位置

名称|值|说明
|---|---|---|
NodeCentered|0|变量存储在网格节点上
FaceCentered|0|变量存储在面心（边界面）
CellCentered|0|变量存储在单元中心

### `FieldIO`

变量输入输出特征

名称|值|说明
|---|---|---|
IncomingDofs|0|变量为输入自由度，如位移、温度
IncomingLoads|1|变量为输入载荷，如力、热流
OutgoingDofs|2|变量为输出自由度，如位移、温度
OutgoingLoads|3|变量为输出载荷，如力、热流

### `FaceTopo`

边界面单元的几何形状

名称|值|说明
|---|---|---|
BAR2|0|2-node 线单元
BAR3|1|3-node 线单元（二次）
TRI3|2|3-node 三角形
TRI6|3|6-node 三角形（二次）
QUAD4|4|4-node 四边形
QUAD8|5|8-node 四边形（二次）
POLYGON|6|N-node 多边形

### `ZoneTopo`

耦合边界（计算域）几何特征

名称|值|说明
|---|---|---|
ZT_POINTS |0|离散点云
ZT_CURVE  |1|曲线边界
ZT_SURFACE|2|曲面边界
ZT_VOLUME |3|实体计算域

### `ZoneShape`

耦合边界（计算域）几何形状

名称|值|说明
|---|---|---|
ZS_POINT    |0|退化为单个点
ZS_COLINEAR |1|所有点共线
ZS_COPLANER |2|所有点共面
ZS_GENERAL  |3|三维曲面或实体

## 2. `VectorInt`

动态数组，元素类型为：`int32_t`

成员函数|原型|功能
|---|---|---|
Vector|Vector()<br>Vector(size)<br>Vector(size,value)|构造函数
.rank||数组的维度，1=向量，2=矩阵，other=多维数组
.numel||数组中元素总数
.empty||数组是否为空（不包含任何元素）
.size||数组尺寸
.clear||删除所有元素
.reserve|reserve(n)|预分配存储空间
.resize|resize(new_size)|更改数组尺寸
.fill|fill(value)|将所有元素值设置为指定值
.swap|swap(other_vec)|与另一个数组交换内存地址和尺寸信息，非常快
.swap_elements|swap_elements(other_vec)|交换两个数组的所有元素值，数组长度必须相同
.copy_elements|copy_elements(src_vec)|将`src_vec`的元素值拷贝到当前数组，数组元素个数必须相同
.push_back|push_back(value)|在数组末尾追加一个元素
.pop_back||删除数组末尾元素


其它用法：
```python
a=VectorInt(3)
a[0]=0  # 支持[]运算符
a[1]=1
a[2]=2
print(a)  # 打印数组内容：{0,1,2}
```

## 3. `VectorIntG`

动态数组，元素类型为：`int64_t`

成语与用法与`VectorInt`相同。

## 4. `Vector`

动态数组，元素类型为：`double`

成员函数和用法同 `VectorInt`，下面给出额外的成员函数：

成员函数|原型|功能
|---|---|---|
.zero||将数组所有元素值置为0
.norm_sq||计算向量模值的平方
.norm||计算向量模值
.min||返回数组最小元素值
.max||返回数组最大元素值
.mean||返回数组元素的平均值
.dot|dot(other_vec)|计算向量内积，数组长度必须相同

## 5. `Matrix`

矩阵，元素类型为：`double`

成员函数|原型|功能
|---|---|---|
Matrix|Matrix()<br>Matrix(m,n)<br>MAtrix(m,n,value)|构造函数
.rank||数组维度，返回2
.numel||矩阵元素总数
.empty||矩阵是否为空
.clear||删除矩阵所有元素
.reserve|reserve(max_size)|预分配存储空间
.resize|resize(m,n,value)|更改矩阵尺寸
.fill|fill(value)|将矩阵所有元素值置为`value`
.swap|swap(other_mat)|与另一个数组交换内存地址和尺寸信息，非常快
.swap_elements|swap_elements(other_mat)|交换两个矩阵的所有元素值，尺寸必须相同，慢
.copy_elements|copy_elements(src_mat)|将`src_mat`的元素值拷贝到当前矩阵，元素个数必须相同
.zero||将矩阵所有元素值置为0
.identity||将矩阵对角元置为1，非对角元置为0
.inverse||将矩阵设置为其逆矩阵，返回值：true=成功，false=矩阵奇异
.apply|apply(x,y)|计算 `y = A . x`，其中`A`为当前矩阵，`x`、`y`为矩阵或向量
.apply_add|apply_add(x,y)|计算 `y = y + A . x`，其中`A`为当前矩阵，`x`、`y`为矩阵或向量

其它用法：
```python
a=Matrix(3,3)
a[0,1]=1 # 将第0行第1列元素设置为1
print(a) # 输出：{{0,1,0},{0,0,0},{0,0,0}}
```

## 6. `Vec3`

三维向量

成员函数|原型|功能
|---|---|---|
Vec3|Vec3()<br>Vec3(x,y,z)|构造函数
.assign|assign(x,y,z)|设置分量值
.fill|fill(value)|填充数值
.swap|swap(other_vec3)|交换两个向量元素值
.x||属性x，用法: a.x=1
.y||属性y
.z||属性z
.norm||计算向量模值
.norm_sq||计算向量模值的平方
.normalize||向量归一化
.dot|dot(other_vec)|计算向量内积
.cross|cross(other_vec)|计算向量外积
.distance|distance(other_vec)|计算与向量`other_vec`的欧式距离
.distance_sq|distance_sq(other_vec)|计算与向量`other_vec`的欧式距离的平方

支持的运算符：`+,+=,-,-=,*float,/float,float*`
代数运算|说明
|---|---|
a+b|计算向量和
a+=b|计算：a=a+b
a-b|计算向量差
a-=b|计算：a=a-b
a*scalar|向量缩放
scalar*a|向量缩放
a*=scalar|向量缩放
a/scalar|向量缩放
a/=scalar|向量缩放

## 7. `IndexSet`

元素编号集合

成员函数|原型|功能
|---|---|---|
IndexSet|IndexSet()<>IndexSet(list)|构造函数
.clear||
.create|create(list)|
.add|add(global_unique_id)|
.contains|contains(global_unique_id)|
.find||
.empry||
.size||
.l2g||
.g2l||
