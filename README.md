# TODO

- [ ] sgRNA和barcode对应
- [ ] 对匹配score进行带温度的softmax，从而把read的count数分配到每个reference
- [ ] mutant和ref只对treat过滤，因为control在merge到treat上时，treat没有的control自然没有。过滤control还会把本来存在的矫正过滤掉。
- [ ] 把kim加回去
- [ ] 仔细调整max_freq_temN
- [ ] 实在不行，保底拟合矫正del长度
- [ ] 利用可视化技术精确定位不好的reads
- [ ] mermaid diagram README for workflow
