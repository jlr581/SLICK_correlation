[t1,y1] = textread("test_series1.txt","%f %f");
[t2,y2] = textread("test_series2.txt","%f %f");
cor_value = correlation(y1,y2,t1,t2)
