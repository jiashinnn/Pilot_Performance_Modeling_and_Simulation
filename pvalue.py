from scipy import stats

# Data
performance_actual = [0.5774, 0.5705, 0.7899, 0.2886, 0.9264,0, 1, 0.7382, 0.4368, 0.7316]
performance_calculated =[0.6498, 0.44626, 0.70471, 0.30208, 0.87893, 0.1438, 0.88786, 0.60491, 0.42483, 0.62837]

# Paired sample t-test
t_statistic, p_value = stats.ttest_rel(performance_actual, performance_calculated)

print("T-statistic: ", t_statistic)
print("P-value: ", p_value)