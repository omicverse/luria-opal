# 暗物质全 Pfam 筛查：结果

## 问题

在去除宿主-Xer 伪影后的 MGE 邻域暗物质（clean_case_v2，79,415 条 orphan）里，用全 Pfam（30,134 模型）`--cut_ga` 筛查是否有未被原始 37 个扫描模型覆盖的蛋白质家族显著富集。

## 设计

- case：clean_case_v2 随机抽 10,000 条（seed=42）；background：bg.faa 随机抽 10,000 条。
- `hmmscan --cut_ga` 对 Pfam-A.hmm（由外部管家用 hmmsearch 完成并转成 hmmscan 口径，本场复核）。
- 每家族 Fisher 富集（case vs bg）+ Benjamini-Hochberg FDR（30,134 次检验）+ 11 管家域中性对照 + 长度分层（>400 桶）。

## falsifier 检查

- **中性对照通过**：11 个管家域全部耗竭（odds 0.156~0.833，中位 0.566），无一 >1.0。
- **长度配对缺陷**：background 未按长度匹配（case 中位 211 aa vs bg 273 aa；case 短蛋白 36% vs bg 23%）。全长度富集偏向短家族；>400 桶是更可辩护的口径。
- **基因组配对缺陷（本轮新发现）**：case 79,415 条 / 72,086 contig，bg 102,336 条 / 97,675 contig，共有 contig 仅 5,757（bg 中 5.9% 来自 case contig）。对比实为「有整合酶的基因组 vs 另一批基因组」，不是「整合酶旁 vs 同基因组别处」。这条缺陷**约束低倍数条目的解释**；高倍数（7–40×）条目在基因组构成 1.2 倍差异下仍站得住（此为判断，未单独用数据文件量化）。

## 结果：32 个 FDR q<0.05 富集家族

**约 28 个是移动元件基因谱**：转座酶（Y2_Tnp ∞、Zn_Tnp_IS91 ∞、rve 39×、DDE_Tnp_Tn3 ∞、DDE_Tnp_IS66 7×、DDE_Tnp_1 4.9×、Transposase_20 9.3×、rve_3 6×）、整合酶配件（Phage_int_SAM_1/3/4/5、Phage_int_M、Zn_ribbon_recom 27.5×）、excisionase（Phage_AlpA 30×）、七个 HTH DNA 结合域、ParB_N 7×、RepSA ∞。

其中 **Phage_integrase 9.88×、Resolvase 27.5×、Recombinase 39.5×、Phage_int_SAM_1、DEDD_Tnp_IS110 本身就是 case 的选择判据**（case 定义为整合酶邻域），是循环论证，不能算发现。

**32 个里另外 3 个不是 MGE 基因谱**：
- **RuvA_N（7.0×，35/5）+ RuvA_C（6.6×，33/5）**：宿主 Holliday junction 修复蛋白（见下）。
- **DUF2924（∞，32/0，q=9.2e-8）**：一个未表征域，但与 RAMA 域（PF18755/IPR040843）部分重叠（见下）。
- **DUF4158（∞，13/0，q=0.021）**：更正——不是未表征孤儿域，是大 Tn3 转座酶的 N 端附属域（13 条中位 984 aa，10/13 与 DDE_Tnp_Tn3 共存），属已计入的 28 个 MGE 基因谱。

Pfam 可命名率：case 63.6%（6,359/10,000），bg 78.4%。完全无 Pfam 命中的残余：case 36.4%、bg 21.6%——这部分是短、进化快的 orphan，正是已报道的「原噬菌体小蛋白」现象。

## 四条开线，全部判清

### 1. RuvA（宿主 ruvAB 操纵子 = 整合热点，非新蛋白）

RuvA_N 7.0×（35/5）、RuvA_C 6.6×（33/5）。33 条同时带两域、长度中位 203 aa（全长宿主型 RuvA，大肠杆菌 RuvA=203 aa）、33 个独立 contig（非共线性伪影）。swissprot 全命中「Holliday junction branch migration complex subunit RuvA」。

**判据（ruvB 邻接）**：33 个 RuvA 位点里 31 个有 ruvB 紧邻（relpos ±1，19 个 +1、12 个 −1），即宿主 ruvAB 操纵子。RuvA 富集是宿主重组修复位点是整合热点，与 pyrG/guaA 同形状，不是 MGE 货舱。

**prior art**：查过但没找到——PubMed 检索 RuvA 作为整合位点，命中全是 RuvA 作为 Holliday junction 解离酶的文献（不同查询式 4~15 命中，无一描述 ruvA/ruvAB 作为整合位点）。最近 SpyCIM1（PMID 26701803）整合位点在 mutL，ruvA 是下游被顺式沉默的旁邻基因，非整合靶。→ 宿主 ruvAB 是整合热点这件事本身未被报道，但它是宿主基因，不是新可编程系统。

### 2. DUF2924（未表征域，32/0）

DUF2924（PF11149）case 32 / bg 0，q=9.2e-8。blastp swissprot 全是弱非特异命中（pident 25–48%、evalue 0.001–4.8，命中核迁移蛋白/coiled-coil/核糖体 uL22B/trigger factor 等无关蛋白），无可名身份。

**prior art（更正）**：上一版按域 ID（`DUF2924`/`PF11149`）查 PubMed 零命中即判「无 prior art」是错判——文献挂在 RAMA 名下。DUF2924 与 RAMA 域（PF18755/IPR040843，clan CL0831/SAND）**部分重叠**（Tan/Aravind/Zhang 2024, Genome Biol Evol 16:evae171, PMID 39106433 明确写「RAMA DNA-binding domain (partly matching Pfam DUF2924)」）；RAMA 由 Iyer/Zhang/Aravind 2016（BioEssays 38:27-40, PMID 26660621）定义为「restriction enzyme adenine methylase associated / predicted modified-DNA reader」。DUF2924 与重组酶/前噬菌体的邻域关联**已发表**：Wolbachia prophage WO 保守排列「large serine recombinase → DUF2924 → ankyrin」（Grève 2024, Front Microbiol 15:1416057, PMID 39238888；Fallon 2023, Insects 14:516, PMID 37367332 亦写「recombinase immediately downstream of DUF2924」）。→ 判定从「无 prior art」改为：**缺少直接的原核分子功能表征，但与预测的 RAMA DNA 结合域部分重叠，且 MGE 关联已被 Wolbachia 文献描述。**

### 3. DUF4158（更正：Tn3 转座酶附属域，非未表征域）

DUF4158（PF13700.12）case 13 / bg 0，q=0.021。**更正**：13 条长度中位 984 aa，仅 3/13 是单域蛋白，10/13 与 DDE_Tnp_Tn3（Tn3 家族转座酶）共存——即它是大型 Tn3 转座酶的 N 端附属域，属于已计入的 28 个 MGE 基因谱，不是独立的暗物质孤儿域。上一版把它与 DUF2924 并列成「两个未表征域」是错的。

### 4. 循环论证家族（非发现）

Phage_integrase / Resolvase / Recombinase / Phage_int_SAM_1 / DEDD_Tnp_IS110 这 5 个是 case 的定义判据，其富集是循环论证，不构成发现。

## 结论：无新可编程系统；一个未表征域（DUF2924，部分匹配 RAMA）是残余线索

**这次暗物质筛查对「新可编程系统」是阴性的**：32 个富集家族里没有一个是新的可编程核酸酶或新已表征的家族。约 28 个是预期的移动元件基因谱（重组/转座机制），RuvA 是宿主整合热点，5 个是循环论证。

**但「阴性」不等于「全部判清」**：DUF2924（32/0）是**一个**未表征域，在整合酶邻域显著富集；它与 RAMA 域（PF18755/IPR040843）部分重叠、与重组酶/前噬菌体的关联已发表（Wolbachia 文献），但原核分子功能尚未直接表征。它是「未知」而非「阴性」——是未来功能表征的候选，但没有任何证据表明它是新的可编程系统。DUF4158（13/0）经更正为 Tn3 转座酶附属域，不再是暗物质信号。

**不配对缺陷写清楚**：本对比的 case 与 bg 几乎不是同一批基因组（共有 contig 5.9%），所以这是一个「有整合酶的基因组 vs 另一批基因组」的对比，不是严格的「整合酶旁 vs 同基因组别处」。高倍数富集（转座酶/整合酶配件 7–40×）不受此缺陷影响而站得住；低倍数条目解释不干净。

## 数据文件

- `family_enrichment.tsv`（4,857 个家族，case/bg/odds/p/FDR/>400 桶/neutral）
- `v2/*.domtbl`（8 块完整 hmmsearch 结果）
- `ruvA_loci_all.faa`（33 contig / 12,876 蛋白，可能被后台进程截断，重建中）、`ruvA_contigs.txt`、`chk_ruvA_full.faa`
- `ruvA_neighbors_swissprot.tsv`（188 邻居 swissprot 命中）
- `mge_machinery_accs.tsv`（142 个整合酶/转座酶/重组酶 Pfam 模型）
