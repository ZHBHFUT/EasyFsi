# EasyFsi python interface

张兵，zhangbing@hfut.edu.cn，2023-06-07

[toc]

## 1. Constants

### 1.1 `FieldLocation`

变量存储位置

名称|值|说明
|---|---|---|
NodeCentered|0|变量存储在网格节点上
FaceCentered|0|变量存储在面心（边界面）
CellCentered|0|变量存储在单元中心

### 1.2 `FieldIO`

变量输入输出特征，见：`Field::iotype`、`FieldInfo::FieldInfo()`

名称|值|说明
|---|---|---|
IncomingDofs |0|变量为输入自由度，如位移、温度
IncomingLoads|1|变量为输入载荷，如力、热流
OutgoingDofs |2|变量为输出自由度，如位移、温度
OutgoingLoads|3|变量为输出载荷，如力、热流

### 1.3 `FaceTopo`

边界面单元的几何形状，见：`Boundary::add_face`、`bd_add_face`、`bd_face_type`

名称|值|说明
|---|---|---|
BAR2   |0|2-node 线单元
BAR3   |1|3-node 线单元（二次）
TRI3   |2|3-node 三角形
TRI6   |3|6-node 三角形（二次）
QUAD4  |4|4-node 四边形
QUAD8  |5|8-node 四边形（二次）
POLYGON|6|N-node 任意多边形（节点数大于4）

### 1.4 `ZoneTopo`

耦合边界（计算域）几何特征，见 `Boundary::topo()`.

名称|值|说明
|---|---|---|
ZT_POINTS |0|离散点云
ZT_CURVE  |1|曲线边界
ZT_SURFACE|2|曲面边界
ZT_VOLUME |3|实体计算域

### 1.5 `ZoneShape`

耦合边界（计算域）几何形状，见 `Boundary::shape()`.

名称|值|说明
|---|---|---|
ZS_POINT    |0|退化为单个点
ZS_COLINEAR |1|所有点共线
ZS_COPLANER |2|所有点共面
ZS_GENERAL  |3|三维曲面或实体

### 1.6 `InterpolationMethod`

插值方法，在 `Interpolator::compute_interp_coeff()` 中使用

名称|值|说明
|---|---|---|
LocalXPS   |0|局部样条
Projection |1|等参逆变换
Mapping    |2|所有点共面
Automatic  |3|三维曲面或实体

## 2. `VectorInt`

动态数组，元素类型为：`int32_t`

成员函数|原型|功能
|---|---|---|
Vector|Vector()<br>Vector(size)<br>Vector(size,value)|构造函数
.rank ||数组的维度，1=向量，2=矩阵，other=多维数组
.numel||数组中元素总数
.empty||数组是否为空（不包含任何元素）
.size ||数组尺寸
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
.zero   ||将数组所有元素值置为0
.norm_sq||计算向量模值的平方
.norm||计算向量模值
.min ||返回数组最小元素值
.max ||返回数组最大元素值
.mean||返回数组元素的平均值
.dot |dot(other_vec)|计算向量内积，数组长度必须相同

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

元素编号集合，见：`Boundary::nodes()`

成员函数|原型|功能
|---|---|---|
IndexSet|IndexSet()<br>IndexSet(list)|构造函数
.clear  |clear()| 清空当前集合
.create|create(list)|从全局编号列表创建集合
.add|add(global_unique_id)|添加一个全局编号
.contains|contains(global_unique_id)|判断是否包含全局编号
.find |find(global_unique_id)|查找全局编号
.empty|empty()|当前集合是否为空
.size |size()|当前集合的元素个数
.l2g  |l2g(local_id)|获取局部编号对应的全局编号
.g2l  |g2l(global_unique_id)|获取全局编号对应的局部编号
[] | `[local_id]`|下标运算符，同`l2g`

## 8. `KDTree`

三维KD-Tree数据结构，见：`Boundary::kdtree()`

成员函数|原型|功能
|---|---|---|
KDTree|KDTree()<br>KDTree(coords,npts,persistent)|构造函数
.clear  |clear()| 清空当前数据结构
.create|create(coords,npts,persistent)|创建KDTree
.search|search(n_pnts,pnts,n_query,idxs,dist_sq,epq)|搜索最近邻居
.empty|empty()|判断当前KDTree是否为空
.size|size()|获取KDTree中的点个数
.data|data()|获取点坐标数组首地址
.bbox_lo|bbox_lo()|获取坐标范围下限
.bbox_ho|bbox_ho()|获取坐标范围上限
.point|point(id)|获取第i个点坐标

## 9. `MeshConnectivity`

CSR格式的网格元素连接信息数据结构，见：`Boundary::face_nodes`.

成员函数|原型|功能
|---|---|---|
MeshConnectivity|MeshConnectivity()|构造函数
.clear  |clear()| 清空当前数据结构
.reserve|reserve(nrow,ndata)| 预分配存储空间，用于加速添加操作
.nrow   |nrow()|获取行数
.ndata  |ndata()<br>ndata(row)|获取（指定行）列元素总数
.empty  |empty()|是否为空
.ia     |ia()|偏移量数组
.ja     |ja()|列编号数组
.push_back|push_back(list)|添加一行
[]      |`[i,j]`|获取第i行第j个列元素

## 10. `Boundary`

组装后的耦合边界数据，见：`Application::add_coupled_boundary()`、`Application::boundary()`

成员函数|原型|功能
|---|---|---|
Boundary|Boundary()|构造函数
.clear  |clear()| 清空边界数据
.reserve|reserve(nnode,nface,nface_nodes)| 预分配存储空间，用于加速添加操作
.add_node|add_node(x,y,z,unique_id)|添加一个节点
.add_face|add_face(type,nodes)|添加一个面
.set_face_cent|set_face_cent(face,fcx,fcy,fcz)|设置面心坐标
.set_face_area|set_face_area(face,farea)|设置面积
.compute_metics|compute_metics()|计算面积、法向、面心的等几何量
.kdtree|kdtree()|获取节点KDTree
.load|load(file)|从文件创建边界
.save|save(file)|保存到文件
.topo|topo()|获取边界拓扑形状
.shape|shape()|获取边界几何形状
.contains_polygon|contains_polygon()|判断是否包含多边形面（五边形以上）
.contains_high_order_face|contains_high_order_face()|判断是否包含二次及以上单元
.all_high_order|all_high_order()|判断是否都是二次及以上单元
.nnode|nnode()|节点个数
.nface|nface()|面个数
.get_field|get_field(name)|获取指定名称的场变量
.node_coords|node_coords(i)|获取第i个节点坐标
.face_centroid|face_centroid(i)|获取第i个面面心坐标
.face_area|face_area(i)|获取第i个面的面积
.face_normal||获取第i个面的单位外法线方向矢量
.face_type||获取第i个面的类型
.register_field||注册一个场变量
.remove_all_field||移除所有场变量

## 11. `Interpolator`

插值器对象。

成员函数|原型|功能
|---|---|---|
Interpolator|Interpolator()|构造函数
.clear|clear()|清空当前插值器
.add_source_boundary|add_source_boundary(bound)|添加源边界（该边界的自由度值是已知的，载荷是未知的）
.add_target_boundary|add_target_boundary(bound)|添加目标边界（该边界的自由度值是待插的，载荷是已知的）
.compute_interp_coeff|compute_interp_coeff(method,max_donor)|计算插值系数，最大邻居应<=20
.save_coefficients|save_coefficients(file)|保存插值系数到文件
.load_coefficients|load_coefficients(file)|从文件读取插值系数
.interp_all_dofs_s2t|interp_all_dofs_s2t()|插值目标边界的所有自由度值
.interp_all_load_t2s|interp_all_load_t2s()|插值源边界的所有载荷值
.interp_dofs_s2t|interp_dofs_s2t(name)|插值目标边界指定名称的自由度值
.interp_load_t2s|interp_load_t2s(name)|插值源边界指定名称的载荷值

## 12. `Communicator`

通信子对象，用于在进程间传递数据。

成员函数|原型|功能
|---|---|---|
.set_constant|set_constant(name,value)|设置参数
.set_function|set_function(name,value)|设置函数指针
.rank|rank()|获取当前进程在该通信子下的编号，从0开始
.size|size()|获取该通信子下的进程个数
.disconnect||断开连接
.send|send(obj,dest,tag)|发送数据到指定进程
.recv|recv(obj,src,tag) |从指定进程接收数据

## 13. `MPICommunicator`

MPI通信子，从`Communicator`派生。

成员函数|原型|功能
|---|---|---|
MPICommunicator|MPICommunicator(mpi_comm,rank,size)|构造函数
.MPI_DATATYPE_NULL||获取或设置常量`MPI_DATATYPE_NULL`
.MPI_INT16||获取或设置常量`MPI_INT16`
.MPI_INT32||获取或设置常量`MPI_INT32`
.MPI_INT64||获取或设置常量`MPI_INT64`
.MPI_FLOAT||获取或设置常量`MPI_FLOAT`
.MPI_DOUBLE||获取或设置常量`MPI_DOUBLE`
.MPI_CHAR||获取或设置常量`MPI_CHAR`

## 14. `SocketCommunicator`

套接字通信子，从`Communicator`派生。

成员函数|原型|功能
|---|---|---|
SocketCommunicator|SocketCommunicator(as_root,np,master_ip,port,timeout_sec)|构造函数


## 15. `Application`

应用程序对象，目前只支持单进程单个对象。

成员函数|原型|功能
|---|---|---|
Application|Application()|构造函数
.clear|clear()|清空当前应用程序对象
.add_coupled_boundary|add_coupled_boundary()|添加一个耦合边界，返回该边界的引用
.boundary_num|boundary_num()|获取耦合边界个数
.boundary|boundary(ib)|获取第ib各边界，编号从0开始
.set_field_function|set_field_function(getter,setter)|设置边界变量读写操作的函数
.register_field|register_field(name,ncomp,location,iotype,units)|注册一个变量
.start_coupling|start_coupling(inter_comm)|开始耦合，使用指定通信子实现不同`Application`之间的数据传递
.exchange_solution|exchange_solution(time,user_data)|插值并交换所有耦合边界物理量
.stop_coupling|stop_coupling()|停止耦合，会待用通信子的`disconnect`函数断开连接
.save_tecplot|save_tecplot(file,without_fields)|保存结果到Tecplot文件

