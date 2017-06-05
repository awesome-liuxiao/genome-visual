# 요약
차세대 심층 시퀀싱 플랫폼에 의해 생성 된 전례없는 양의 데이터와 대규모 데이터 분석에 대한 수요가 증가함에 따라 소프트웨어 응용 프로그램을 최적화하는 것이 필수적입니다. 따라서 정렬 정확도와 계산 비용 사이의 균형이 가장 적합한 MSA 프로그램의 중요한 지표가 되었습니다. 우리는 벤치 마크 정렬 데이터 세트 [BAliBASE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC29792/)에 대해 `CLUSTALW`, `CLUSTAL OMEGA`, `MUSCLE`의 3 가지 MSA 프로그램의 정확성과 비용을 비교하고 포함 된 일부 구현의 관련성에 대해 논의했습니다. 각 프로그램의 알고리즘에서 정렬의 정확도는 BAliBASE에서 제공 한 두 가지 표준 채점 함수, 쌍의 합 및 전체 열 점수를 사용하여 계산되었으며 계산 비용은 최고 메모리 사용량과 실행 시간을 수집하여 결정되었습니다.

# 도입

서열정렬(sequence alignment)은 단백질서열(protein sequence)이나 핵산서열(nucleic acid sequence) 사이의 상관관계를 나타내는 것으로 서열들이 기능적으로 혹은 진화적으로 어느 정도 연관성이 있고, 서열의 어느 부분들이 그러한 연관성을 가지고 있는가를 나타내는 방법으로 표시될 수 있다. 서열을 정렬한 후에 분석하는 것은 두 생물을 비교할 때 머리-머리, 다리-다리, 꼬리-꼬리와 같이 공통적인 기원을 갖는 기관을 비교대상으로 선택하는 것과 같은 의미를 포함한다. 즉, 서열정렬은 어떤 기준에 의한 계산으로부터 만들어지는데 그 계산 과정이 복잡하고 계산량이 많기 때문에 컴퓨터 프로그램을 사용해야만 한다.

서열정렬은 1) 관심 대상인 sequence와 상동성(homology)이 높은 서열들을 알아내어 서열의 기능을 유추하거나, 2) 관련 있는 서열들간의 정량적인 상관관계(진화적인 연관성과 같은)나 관련 기능 부위 등을 예측하기 위한 목적으로 이용된다. 또한 3) 특정 집단에 적용할 유전자탐침(DNA probe) 또는 핵산증폭반응을 위한 universal primer의 제작 등에도 서열정렬이 필요하며, 4) 유전자의 전체염기서열을 얻기 위한 contig(서열 조각) 구성작업에도 서열의 정렬이 필요하다. 

# 관련 연구
## 상동성(homology) 분석의 기초

서열비교(sequence comparison)는 서열간의 관계를 나타내기 위해 비교하는 것인데, 이 때 서열이란 단백질서열이나 핵산서열을 말한다. '서열간의 관계'는, 서열들이 기능적으로나 진화적으로 어느 정도 연관성이 있고, 서열의 어느 부분들이 그러한 연관성을 가지고 있는가를 나타내는 것이며, 이것을 표현하는 것을 서열정렬(sequence alignment)이라고 한다. 

서열간의 기능적이거나 진화적인 연관성을 상동성(homology)이라고 하는데, 생물학적인 그 의미를 서열분석에서 명확히 정의하기는 어렵고, 경우에 따라 조금씩 그 정의를 달리한다. 서열비교는 서열정렬로 표시하므로, 일반적으로 서열비교는 서열정렬과 같은 의미로 사용한다. 
```
>seq1 MAVDPPKADPKGSSGGFHRKSADDPKNGGPALGSREDQSAKAGCCDPPKADPKGGSRDRVRRCIRANLLVASLLTV
>seq2 DPPKADPKGGGFHRNARKKDNGGPALGSREDQSAKAGGCCQWERTASGSRDRVRRCIRANLLVLLTVAAVV

pair alignment
> seq1 MAVDPPKADPKGSSGGFHR ...
          |||||||||  |||||
> seq2 ---DPPKADPKG--GGFHR ...
```

이러한 정렬은, 어떤 기준에 의한 계산으로부터 만들어지고 그 계산 과정이 복잡하고 계산 량이 많기 때문에 컴퓨터 프로그램을 사용해야만 한다. 사실상 모든 서열정렬은 프로그램에 의해 이루어지므로, 서열정렬을 하여 서열을 분석하고자 하는 생물학자는 프로그램의 사용법과 해석 방법을 습득해야 한다.

## 서열정열의 종류

서열의 갯수
- 쌍정렬(pairwise alignment) - 2개의 서열에 대한 정렬
- 다중정렬(multialignment) - 3개 이상의 서열에 대한 정렬

상동성의 종류
- 전역정렬(global alignment) - 전역상동성(global homology)에 대한 정렬
- 지역정렬(local alignment) - 국부상동성(local homology)에 대한 정렬

위에서 예로 보인 정렬은 쌍정렬에 해당하고, 정렬에서 사용된 서열이 3개 이상인 경우를 다중서열정렬(multiple sequence alignment), 혹은 간단히 다중정렬(multialignment)이라고 한다. 2개의 서열에 대한 쌍정렬도 다중서열정렬의 한 경우로 포함할 수 있고 실제로 다중서열정렬 프로그램을 사용하여 쌍정렬을 생성할 수도 있지만, 다중서열정렬의 생성 방법이 쌍정렬에서 사용한 계산과정보다 훨씬 복잡한 단계를 거치고 계산방식이 다르기 때문에 3개 이상의 서열로 구분하여 이야기하는 것이 일반적이다. 다중서열비교(multiple sequence comparison)는 역시 다중서열정렬을 사용하여 표현되며, 쌍서열비교(pairwise sequence comparison)와 그 응용분야나 분석의 목적이 다르다. 

쌍정렬(pairwise alignment) - 쌍정렬(pairwise alignment)은 정렬대상이 2개의 핵산-핵산 또는 단백질-단백질 서열이며, 주로 상동성(homology)을 나타내기 위해 사용된다. 상동성 검색을 위한 가장 기초적인 방법이므로, 보다 복잡한 문제의 알고리즘들은 대부분 쌍정렬 알고리즘과 연관되어 있다. 

다중정렬(multialignment) - 다중정렬(multialignment)은 다중서열정렬(multiple sequence alignment)이라고도 하며 3개 이상의 핵산 또는 단백질 서열을 하나의 정렬로 나타내는 것으로, 패밀리(family) 분석, 계통관계(phylogeny) 분석, 도메인(domain) 분석 등의 기능분석(functional analysis) 연구를 위해 매우 다양하게 사용된다. DNA 염기서열 결정이 급속히 진전되면서 같은 그룹의 여러 개의 핵산이나 단백질 서열에 대한 비교가 필요하게 되고, 단백질 패밀리나 관련 핵산의 모티프(motif) 서열 등에 대한 자료가 증가하고 기능 분석 연구가 활발해지고 있으므로, 효율적인 다중정렬의 필요성이 지속적으로 제기되고 있다. 

전역정렬(global alignment) - 서로 동일한 종류의 핵산이나 단백질 sequence이어서 비교하는 sequence를 전체적으로 비교하여 최대의 상동성이 나타나도록 alignment하는 경우, 그러한 상동성을 전역상동성(global homology)이라 하며, 그 정렬을 전역정렬(global alignment)이라고 한다. 동일한 기능을 가진 단백질이나 핵산으로 높은 global homology를 나타내는 경우에는 local homology도 global homology를 보이는 부분과 유사한 부분을 나타내게 된다. 그러나 전역정렬과 국부정렬은 서로 다른 기준에 의해 다른 방식의 계산을 해서 다른 결과를 보여주기 때문에, 적용된 알고리즘이 서로 다르고, 사용되는 프로그램도 구별된다. 

지역정렬(local alignment) - 두 서열의 어떤 부분 서열이 높은 상동성을 가지는가를 나타내기 위해 정렬하는 경우, 그러한 상동성을 국부상동성(local homology)이라 하고, 그 정렬을 국부정렬(local alignment)이라고 한다. 기능적으로 일부분만이 관련되어 있으며, 진화적으로도 일부분만이 상동성을 보이는 경우에는, local alignment를 실행하는 것이 효과적일 것이다. 또한 유사한 기능의 sequence이나 global homology가 낮아서 쉽게 상동성을 나타내지 못하지만, 일부분에서 높은 local homology가 있는 경우도 local alignment를 실행하는 것이 효과적이다. 따라서 많은 경우 local alignment가 상동성 sequence를 검색하기 위한 용도로 사용된다. 

## 현재 연구

1) 쌍정렬(pairwise alignment)

1970년 Needleman과 Wunsch의 dynamic programming 기법에 의한 pairwise alignment algorithm이 발표된 이후, Smith와 Waterman 등이 일반적인 길이의 gap에 대해 이 algorithm을 확장함으로써, 두 sequence 간의 alignment를 위한 실용 가능한 프로그램들이 구현되기 시작하였다. 전산학의 graph algorithm 이론이나, pattern matching algorithm 등이 주로 생물학 패턴의 검색에 많이 응용되어 여러 응용 알고리즘들이 개발되어 왔으며, 데이터베이스의 개발과 함께 대량의 데이터를 대상으로 하여 빠른 실행 시간이 요구되어 여러 heuristic algorithm들이 개발되었다.

Lipman, Wilbur, Pearson 등이 hashing과 window 설정의 방법을 사용하여 고안한 FASTP/FASTN는 local homology alignment를 포함한 FASTA로 발전하여 데이터베이스 대상의 homology search에 사용되었고, Altschul과 Karlin이 homologous candidate fragment 선정과 extension 방법을 사용하여 개발한 BLAST는 데이터베이스 대상의 local homology search를 위해 가장 많이 활용되어 왔으며, 최근에는 BLAST의 결정적인 단점으로 gap 처리가 안되던 것을 해결한 BLAST2가 발표되어 BLAST를 대체하고 있다. 같은 목적으로 MPP 기종에서 병렬 프로그래밍된 MPsrch가 상용 프로그램으로 개발되어 있다. 

2) 다중정렬(multialignment)

1985년 Murata 등은 서열쌍정렬 방법을 3개의 서열에 확장하여 적용시킨 THREE를 발표하였으나 실행시간과 메모리 필요량이 O(LN)(길이 L인 N 개의 서열에 대해)이어서 실용적이지 못하였다. 즉, N개의 서열에 대해 이 방법을 그대로 확장할 경우, 각각 LN에 비례하는 실행시간과 메모리가 필요하게 되어 사실상 실행이 불가능한 것이다. 따라서, 이를 극복하기 위한 heuristic 알고리즘들이 발표되어 왔다. 1986년 Sobel 등은 substring matching에 의한 pairwise alignment 방법을 확장한 알고리즘을 발표하였고, Bains는 template sequence와의 pairwise alignment의 결과를 dynamic하게 반영하는 알고리즘을 개발하였고, Taylor는 3차구조가 밝혀진 단백질에 대해 같은 그룹에 속하는 단백질의 sequence를 구조 정보를 이용하여 align하는 알고리즘을 발표하였다. Barton이 pairwise alignment에 근거해 결정한 order에 따라 한 sequence씩 merge해 나가는 방법을 발표한 이후 ordering과 merging 방법의 다양한 조합에 의한 algorithm들이 계속 발표되었다.

실용적인 프로그램으로서 gap 처리와 ordering에 있어서 비교적 안정적인 결과를 보인 ClustalV가 Higgins 등에 의해 발표되었고, gap cost를 dynamic하게 처리하는 방법을 도입하여 개선한 ClustalW가 발표되어 발표되었다. 최근에는 사용자 인터페이스를 개선한 ClustalX도 사용되고 있다.

# MSA 정렬방법

Alignment 방법의 일반적인 원칙은 두 쌍의 염기서열의 유사도가 가장 크게 되도록 배열하는 것이다. 즉 mismatch와 gap이 최소가 되도록 한다. 배열 중에 필요하면 삽입(insertion)과 결손(deletion)을 가정하는 공백(gap)을 추가할 수 있으나, 일반적으로 gap penalty를 mismatch 보다 크게 한다.

정렬을 위한 상동성 계산에서는 scoring matrix와 gap cost 등과 같은 몇 가지 parameter를 공통적으로 필요로 한다. 현재의 alignment 알고리즘과 프로그램들의 지능화 수준을 고려할 때, 이들 파라메터의 의미를 파악하고 이에 맞추어 적용하는 것이 보다 의미 있는 alignment를 구하기 위한 시작인 것이다. 

최근까지 발표된 대부분의 다중정렬 알고리즘들은 공통된 틀을 가지고 있다. 대부분의 다중정렬 알고리즘은 서열쌍정렬 알고리즘에서 사용된 기법을 반복하는 방식을 취하며 다음과 같이 구성된다. 

[1] 각 서열쌍에 대해 서열쌍정렬 알고리즘을 사용하여 상동성 값을 구한다.
[2] ①의 결과로부터 다중정렬해 나갈 순서를 정한다.
[3] ②에서 결정된 순서와 정렬 방법에 의해 하나의 공통서열이 형성될 때까지 정렬해 나간다. 

[1]의 과정은 서열쌍정렬의 몇 가지 알고리즘 중에서 하나가 선택적으로 사용된다.
[2]에서 정렬 순서는 각 알고리즘에 따라 다양하게 결정된다. 상동성이 우선적으로 높은 것을 먼저 선정하여 클러스터링(clustering;묶음 지어 나가는 방법)하는데, 클러스터(cluster;묶음 지어진 개체)들간의 상동성을 비교하는 방법이 다양하기 때문에 최종적으로 클러스터링되는 결과가 다양하게 나타날 수 있다.
[3]에서 정렬 방법은, 다음 클러스터링의 정렬에서 기존 클러스터내의 정렬된 서열을 프로파일(profile)로 간주하는 경우와 하나의 서열만을 결정하여 정렬하는 경우 등과 같이 여러 가지 방법들이 사용된다. 

## Scoring matrix

서열정열에는 두 sequence 간의 상동성 계산이 필요하며, alignment에서 연결된 모든 염기쌍, 혹은 아미노산쌍의 상동성의 합으로 계산되는데, 각 염기쌍과 아미노산쌍에 대해 상동성 값을 미리 정해둔 것이 바로 scoring matrix이다. 따라서, 각 alignment를 구하는 과정에서 scoring matrix는 alignment의 형태와 선택에 있어 결정적인 역할을 하게 된다.

염기쌍: 핵산 서열내의 염기간의 치환에 대한 통계는 일반적으로 사용되기 어려운 점이 있기 때문에, 염기쌍에 대해서는 같은 염기쌍 여부만을 1, 0 등으로 표시하는 간단한 방식을 사용한다. 그러나 purine, pyrimidine 처럼 IUB 방식으로 표기되는 불완전한 정보의 비교에 score matrix를 사용하기도 한다.

아미노산쌍: 아미노산쌍에 대해서는 생물학적인 모델과 통계적인 방법에 근거하여 유도한 값을 사용한다. 서열정렬 프로그램에서 옵션(option)으로 선택하는 scoring matrix는 일반적으로 단백질 서열을 비교정렬하는 경우에 필요한 아미노산쌍에 대한 scoring matrix 이며 PAM, BLOSUM이 사용된다. 

### Scoring Matrix 관련 연구

1969년 M. Dayhoff가 mutation data matrix인 PAM을 발표하였고, McLachlan은 아미노산의 genetic, structural similarity 등을 분석하고 alternative amino acids 방법에 근거한 AAAM matrix를 개발하였다. Altschul은 Information theory에 근거한 matrix를 발표하였다. 그러나 실용적인 면에서 PAM이 alignment에서 좋은 결과를 보여주었다. 1992년 S. Henikoff와 J. Henikoff가 protein block으로부터 구성한 matrix인 BLOSUM이라는 새로운 matrix를 발표하고, 그 효용성이 인정되면서, 최근에는 PAM과 함께 표준적인 scoring matrix로 사용되고 있으며, 여러 alignment 관련 프로그램에서 PAM에 우선하여 표준적인 scoring matrix로 채택되고 있다. 

### Scoring Matrix 선택

하나의 scoring matrix를 선택하는 것은, 대부분의 프로그램에서 하나의 matrix를 사용하기 때문에 생긴 제약일 뿐이고, '일반적인 혹은 표준적인 스코링 매트릭스'라는 것은 별로 의미가 없다. 따라서, 사용자는 허용되는 여러 종류의 스코링 매트릭스를 사용하면서, 어느 경우에 더 만족할 만한 정렬이 나타나는 가를 비교하여, 최종적으로 하나의 정렬을 선택하여야 한다.

PAM(Accepted Point Mutation) matrix는 일정한 진화시간 동안의 아미노산간의 치환정도로부터 계산된다. PAM 시리즈는 PAM20, PAM250 등과 같이 PAM(N)으로 표현한 여러 매트릭스인데, 이때 N은, PAM을 유도할 때 사용된 '단위 진화거리'의 배수를 나타낸다. 즉, PAM200은 PAM100의 진화거리의 두 배의 진화거리에 대해(즉, 진화시간에 대해), 각 아미노산쌍간의 연관성이나 상호교환정도를 나타낸다.

Alignment 프로그램에서 PAM을 선택할 경우, 비교하는 두 sequence의 진화시간이 큰 경우에는 높은 N값의 PAM을, 진화적 거리가 작은 경우에는 작은 N값의 PAM을 선택하는 것이 비교적 정확한 alignment를 기대할 수 있다. 한편, 실제로 기능과 관련 있는 단백질 block들로부터 아미노산의 치환정도를 계산하여 구성된 BLOSUM 시리즈도 BLOSUM62, BLOSUM90 등과 같이 BLOSUM(N)으로 표시되는 여러 매트릭스이다. 그러나 BLOSUM(N)은 진화의 단위대신에 N % 이상의 상동성을 가진 sequence segment로부터 clustering된 block 들로부터 계산되었다는 의미이다. 따라서 BLOSUM에서는 N이 클수록 상동성이 큰 서열쌍의 정렬에 적합하다고 할 수 있다. 

## Gap cost

정렬결과에서 '-'로 표시되는 갭(gap)은 진화상 삽입(insertion) 또는 결손(deletion)에 해당하는데, 한 서열에서 생략되어 있거나 다른 서열에서 추가되어 있는 부분을 나타내는 것이다. 정렬 알고리즘에서는 scoring matrix로 아미노산쌍이나 염기쌍에 대해 일정한 값을 주듯이, 이들 gap이 발생하는 경우에 대해서도 일정한 값을 주어, 상동성을 계산하고 적절한 형태의 정렬이 구해지도록 한다. Alignment에서 gap 발생에 대한 cost(penalty)를 얼마의 값으로 두느냐에 따라 gap 발생 장소와 길이, 횟수가 결정되어 alignment의 형태에 역시 중요한 영향을 주게 된다. 

[1] Gap cost 계산모델

Alignment에서 gap 발생의 효과는 gap penalty로 표시하는데, 여러 가지 모델들이 제시되고 구현되었다. 가장 간단한 것이 길이에 상관없이 모든 gap에 대해 상수인 GP를 설정하는 것이다.

Model ① GP = Constant

이에 비해 생물학적인 의미를 좀 더 부여해서, 길이가 긴 gap일수록 짧은 gap에 비해 insertion, deletion의 사건이 반복되었을 가능성이 큰 것으로 가정하여, 길이에 따라 증가하는 모델들이 제시되었다.

Model ② GP=g*L (L은 gap의 길이, g는 상수)

Model ③ GP=GOP+GEP*L (L은 gap의 길이, GOP는 gap opening penalty, GEP는 gap extension penalty)

Model ④ GP(L)=GP(L-1)+gL (gL< gL-1)

가장 일반적으로 사용되는 모델은 Gotoh가 제안한 affine gap cost라 불리는 GP=GOP+GEP*L 의 모델이다. Gotoh의 알고리즘은 계산량의 복잡도가 증가하지 않고 기존의 계산 방식에 약간의 변경을 준 것으로 매우 실용적인 알고리즘으로 평가되고 있다. 

[2] Gap cost의 선택

Gap비용을 객관적으로 설정할 기준은 현재 없다. Gap의 생물학적 의미는 하나의 sequence에서 insertion이 발생했거나, 혹은 다른 sequence에서 deletion이 발생했다는 것이다. 실험실에서의 경험적인 통계와 계통관계 연구로부터 이러한 insertion/deletion은 다른 아미노산(핵산에서는 염기)로의 치환보다는 드문 사건이라고 생각할 수 있다. 따라서, 이것은 치환의 경우보다 상대적으로 낮은 값, 즉 전체 상동성값을 감소시키는 비교적 큰 절대치의 음수값을 준다. 엄격히 말하면, insertion과 deletion에 대해서는 다른 값을 주어야 하지만, 두 sequence를 비교할 경우, insertion과 deletion을 판정할 수가 없으며, 판정이 가능하다 해도(multialignment의 경우는 예측은 할 수 있다), 두 값에 대해 서로 상대적인 값을 줄 근거(계산방법)가 현재 불충분하다.

GOP는 GEP보다 큰 것이 일반적이고(둘 다 음수이므로 절대치로서 크다), GOP와 GEP는 스코링 매트릭스에서의 아미노산쌍 값 중 작은 값보다도 비교적 더 작은 값(음수 절대치로서는 큰 값)을 가져야 하는 것이 일반적이다. 갭 모델은, 스코링 매트릭스에 비해서도 생물학적인 의미를 잘 반영하지 못하고 있다. 그런 면에서 특히, 갭비용의 경우는 서열정렬 프로그램 실행시에 허용되는 범위 내에서 여러 가지 값을 설정하여 실행하는 것이 필요한 것이다. 

==== TODO ====

## ClustalW

## ClustalOmega

## MUSCLE

==== Ref ====
Assessing the efficiency of multiple sequence alignment programs

ClustalW
J. D. Thompson, D. G. Higgins, and T. J. Gibson (1994).CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. Nucleic Acids Res., 22(22):4673 4680. DOI: 10.1093/nar/22.22.4673.

ClustalOmega
F. Sievers, A. Wilm, D. Dineen, T. J. Gibson, K. Karplus, W. Li, R. Lopez, H. McWilliam, M. Remmert, J. Söding, J. D. Thompson, and D. G. Higgins (2011) Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol. Syst. Biol., 7:539. DOI: 10.1038/msb.2011.75.

MUSCLE
R. C. Edgar (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5(5):113. DOI: 10.1186/1471-2105-5-113.
R. C. Edgar (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res., 32(5):1792 1797. DOI: 10.1093/nar/gkh340.
