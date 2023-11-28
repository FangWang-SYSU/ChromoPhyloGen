#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File   : estimateDER.py
@Author : XinWang
"""

import pandas as pd
import statsmodels.api as sm


def estimate_der(tree_obj, cell_relation):
    pseudo_res = cell_relation.__copy__()
    pseudo_res.set_index(pseudo_res['name'], inplace=True)
    brother_node = tree_obj.edges.__copy__()
    brother_node.set_index(brother_node['c'], inplace=True)
    brother_node.loc['root',:] = ['root', 'root', 0]
    brother_node = brother_node.loc[pseudo_res['name'], :]
    pseudo_res[['brother_dd_loc', 'brother_ad_loc']] = pseudo_res.loc[brother_node['p'], ['Mitosis_dd_loc', 'Mitosis_ad_loc']].values

    model_dataset = pseudo_res.dropna()

    train_data = pd.DataFrame({
        'parent_loc': model_dataset['Parent_gain_loc'].astype(float) + model_dataset['Parent_loss_loc'].astype(float),
        'parent_cn': model_dataset['Parent_gain_cn'].astype(float) + model_dataset['Parent_loss_cn'].astype(float),
        'brother_loc': model_dataset['brother_dd_loc'].astype(float) + model_dataset['brother_ad_loc'].astype(float),
        'pseudotime': model_dataset['Pseudotime_tree'].astype(float)
    })

    train_label = model_dataset[['Mitosis_dd_loc', 'Mitosis_ad_loc', 'Mitosis_time']].astype(float)
    train_label = train_label.sum(axis=1)

    train_data['label'] = train_label
    train_data = train_data.fillna(0)

    model_dataset = pseudo_res.loc[~pseudo_res.index.isin(model_dataset.index),:]
    model_dataset = model_dataset.loc[~model_dataset['Root_gain_loc'].isna(), ]
    test_data = pd.DataFrame({
        'parent_loc': model_dataset['Parent_gain_loc'].astype(float) + model_dataset['Parent_loss_loc'].astype(float),
        'parent_cn': model_dataset['Parent_gain_cn'].astype(float) + model_dataset['Parent_loss_cn'].astype(float),
        'brother_loc': model_dataset['brother_dd_loc'].astype(float) + model_dataset['brother_ad_loc'].astype(float),
        'pseudotime': model_dataset['Pseudotime_tree'].astype(float)
    })

    Y = train_data[['label']]
    Y['label2'] = tree_obj.node_data.shape[1] - train_data[['label']]
    glm_model = sm.GLM(
        Y,
        sm.add_constant(train_data[['parent_loc', 'brother_loc', 'pseudotime']]),
        family=sm.families.Binomial()
    )
    glm_model = glm_model.fit()
    glm_pred = glm_model.predict(sm.add_constant(test_data[['parent_loc', 'brother_loc', 'pseudotime']]))
    return glm_pred
