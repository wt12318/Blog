---
title: 【Ch1】 Hands on Machine Learning
author: wutao
date: '2021-01-13'
slug: hands_on_ML
categories:
  - reading notes
tags:
  - ML
  - notes
image : "image1.png"
---

第一章：The Machine Learning Landscape
对机器学习的概述，包括机器学习系统的分类，机器学习的挑战和一些术语的解释(验证集，训练集，测试集，正则化，超参数，交叉验证等)


机器学习对数据挖掘的作用除了预测之外：我们也可以通过检查机器学习模型从数据中学习到的模式来对问题有更深入的理解，也就是机器学习可以帮助人去学习

机器学习系统的类型：

- **监督学习和非监督学习** 在监督学习中有两类算法容易混肴：Anomaly detection and novelty detection 异常检测和新奇检测 异常检测指的是来了新数据，检测其是否属于正常的数据(训练集的数据)，而新奇检测则是检测其是否不同于训练集中的数据 比如训练集中有几千张狗的图像只有1%的图像是吉娃娃，那么来了一个新的吉娃娃图像，那么就判断其为异常数据 而新奇检测则不会将其判断为是新奇的数据

- **Batch and Online Learning**  

  - *batch learning* 不能够增量学习 必须使用所有的数据进行训练 每次有新的数据来的时候必须要和旧数据进行整合来重新训练模型 也就是说当新数据来的时候我们必须要将旧的模型关闭下线，重新训练后再上线运行，所以这种模式也叫做离线学习 offline learning

  - *online learning* 可以依次传入数据进行**增量学习**，传入的形式可以是单个的数据也可以是小批量数据(mini batch) 因此有两种应用场景：依次传入的数据；非常大的数据(拆分进行训练)

    ![image-20210111183408092](https://picgo-wutao.oss-cn-shanghai.aliyuncs.com/img/image-20210111183408092.png)

    ![image-20210111183443429](https://picgo-wutao.oss-cn-shanghai.aliyuncs.com/img/image-20210111183443429.png)

- **Instance** **Based Versus Model** **Based Learning**  依据泛化方式的不同进行分类

  - 基于实例的：根据相似性度量
  - 基于模型：依据训练数据选择模型---训练模型---进行推断(在新的数据上应用，也就是泛化)

机器学习的主要挑战：

- 训练数据的数量不足

- 训练数据的代表性不足 可能有两个原因导致：训练数据比较小(*sampling noise* :由于随机性产生的不具代表性的数据)；尽管数据比较大，但是由于抽样方法的错误导致的抽样偏差(*sampling bias*)

- 低质量数据

- 不相关的特征 *garbage in, garbage out* 机器学习的重要步骤就是我们要选取一系列好的特征来训练我们的模型，这个步骤也叫做**特征工程**( *feature engineering*),包括以下几个步骤：

  - 特征选择 在已有的特征中选择最有用的特征来训练模型
  - 特征提取 组合已有的特征产生更有用的特征，比如使用降维算法
  - 通过收集新的数据来创造(发现)新的特征

- 在训练数据上的过拟合 过拟合指的是模型在训练集上表现比较好，但是在测试集上表现不好，也就是泛化能力不足；解决过拟合可能有以下几种方法：

  - 简化模型：选择更少的模型参数，减少训练数据的属性数量(特征)，约束模型(正则化)
  - 收集更多的训练数据
  - 减少训练数据的噪音：消除数据错误，移除数据中的离群点等

  这里比较重要的模型的约束，也就是**正则化**(*regularization*)，比如在简单线性回归中我们可以限制直线斜率的大小的变化从而限制模型的自由度；在学习过程中应用的正则化的数量是通过**超参数**(*hyperparameter*)来控制的,*超参数是学习算法的参数而不是模型的参数*，需要在训练模型之前就确定好

- 在训练数据上的欠拟合(underfitting) 当模型过于简单而不能学习到数据的内在结构时，那么该模型即使在训练集上也表现不好，可能通过以下方法来改进：

  - 选择一个更复杂的模型，使用更多的参数
  - 通过特征工程喂给学习算法更好的特征
  - 减小模型的约束，比如减少正则化超参数等

评估模型的泛化能力，即模型在没有见过的数据上的表现. 可以使用将数据分为**训练集(*training set* )和测试集(*test set*)** 并在训练集上训练模型，在测试集上测试模型. 模型在新的样本上的错误率叫做泛化误差(generalization error或者out-of-sample error),通过在测试集上评估模型可以估计这个错误率，从而告诉我们模型在没有见过的数据上表现如何；如果模型在测试集上误差很小，而泛化误差比较大，说明该模型在训练数据上过拟合 【一般将数据的80%作为训练集，20%作为测试集，不过这也看数据集的大小而定】

当我们在不同的模型之间做选择的时候，可以看这些模型在测试集上的误差从而选择表现最好的模型

当我们需要在一个模型上进行超参数的调试的时候，比如有100个需要尝试的超参数的值，我们用这100个不同的超参数训练模型并在测试集上评估，选择一个表现最好的超参数，但是这可能会发生一种情况：这个模型在测试集上泛化误差比较低，但是应用到新数据上表现又不是很好 我的理解是：在测试集上反复评估模型调整超参数，从而产生一个在测试集上表现最好的模型，所以这个时候测试集变成了“训练集”，这个模型可能形成了新的过拟合现象

一个普遍的解决方法是从训练集中拿出一部分数据来评估备选的模型并选择最好的模型，这一部分数据叫做***validation set*** 也就是说流程为：在剩下的训练集来训练不同的超参数构成的模型，然后选择在*validation set*上表现最好的模型，这就相当于在全部的训练数据上来得到最好的模型，最后使用测试集得到泛化误差来评估这个模型

还有一个问题就是验证集的大小，如果验证集太小，模型评估就会不精确；如果验证集太大，那么剩下的训练集就会比较小(也就是用来训练模型的数据比较少，一个比喻就是选择短跑健将来参加马拉松). 解决这种问题可以使用**交叉验证(cross-validation)** :使用很多小的验证集，每个模型在余下的训练集上训练后再在验证集上进行评估，对每个验证集都进行这个步骤，最后将一个模型在所有验证集上的评估的平均作为该模型表现的更精确的衡量




