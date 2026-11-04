# 공개 DNA 서열 / 체세포 변이 자료 실사

기준일: 2026-11-04. **메틸화 자료 없이 A/C/G/T 서열 자체를 얻거나 종양의 변이를 확인하려는 경우**의 자료 출처.
원시 WGS, 처리된 변이 호출, 표준 참조 서열은 서로 다른 자료이며 같은 이름으로 취급하면 안 됨.
등록정보와 제공 설명은 확인했지만, 본 패키지에서는 대용량 FASTQ/BAM 다운로드 및 실제 사용자 권한/전송 성공을 검사하지 않았음.

## A. 유전자 이름만 입력해서 바로 받는 참조 서열

| 제공처 | 무엇을 얻나 | 접근 | 주의 |
|---|---|---|---|
| Ensembl REST | 지정 종·유전체 좌표의 reference DNA | 공개 API | 환자 암세포 DNA 아님. 한 요청 최대 10 Mb로 문서화됨 |
| UCSC Genome Browser API | `hg38`, `hg19` 등 reference assembly의 구간 서열 | 공개 API | 0-based half-open 좌표를 사용하는 API 호출에 주의 |

공식 API: https://rest.ensembl.org/documentation/info/sequence_region
https://genome.ucsc.edu/goldenPath/help/api.html

이들은 앱의 종별 유전자 CGR에 알맞지만, 암/정상 **환자별 서열 비교** 자료는 아님.

## B. 실제 공개 종양·정상 WGS 원자료

| 자료 | 유형 | 구체적 공개 입구 | 분석 대상 | 규모/한계 |
|---|---|---|---|---|
| NIST GIAB HG008 | pancreatic cancer tumor/normal WGS + somatic benchmarks | https://www.nist.gov/programs-projects/cancer-genome-bottle | 암·정상 구간 서열 재구성, 분석 검증 | 표준 검체(대규모 임상 코호트 아님); 공개 FTP 데이터 |
| NIST GIAB HG009 | pancreatic cancer liver metastasis-related reference WGS | 같은 NIST 페이지 | 추가 암 기준 검체 | 정상·종양 대응/사용 가능한 변이 벤치마크 파일을 분석 전에 확인 필요 |
| SEQC-II HCC1395/HCC1395BL | breast cancer / matched normal cell lines, WGS FASTQ and variant call sets | NCBI SRA study `SRP162370` · https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/ | 암 변이 반영 CGR을 만들기 좋은 첫 실험 | 환자 임상 코호트가 아니라 cell line benchmark |
| PacBio Revio HCC1395/HCC1395-BL | same tumor/normal cell lines, HiFi WGS reads + methylation calls | https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/ | 실제 DNA 경로와 메틸화 후속 통합 | 자료 용량 매우 큼; 같은 세포주의 다른 실험 |
| CASTLE HCC1395/HCC1395BL | same tumor/normal cell lines, ONT long reads | https://github.com/CASTLE-Panel/castle | long-read 서열/변이 탐색 | SRA 개별 접근번호 `SRR28305167`, `SRR28305166`; 임상 코호트 아님 |

위 데이터는 **원자료 자체가 공개된 사례**를 위주로 작성함. 가용 용량, FTP/SRA 전송 제한, 시퀀싱 검체의 대응 관계는 실제 분석 전에 재확인해야 함.

## C. 넓은 암종 범위의 공개 *변이 결과* (원시 WGS가 아님)

- **NCI GDC open somatic MAF**: TCGA 여러 암종에 걸친 체세포 변이 호출표를 다운로드할 수 있음. 공개 Somatic MAF는 민감한 germline 변이를 제거한 자료로, WES 기반 사례도 포함됨. **WGS 전체 염기서열이라고 표기하면 틀림.** https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format/
- **PCAWG 일부 공개 consensus variant 결과**: 전장유전체 분석에서 추린 체세포 변이·서명 결과. 환자별 원시 WGS FASTQ/BAM과는 별개이며, 일부 자료는 접근 제한임. 기존 ICGC Data Portal은 **2024년 6월 종료**되었으므로 오래된 `dcc.icgc.org` 다운로드 URL을 즉시 유효하다고 가정하면 안 됨. https://docs.icgc.org/ · https://www.nature.com/articles/s41467-020-16785-6

## 우선순위

1. **참조 CGR 만들기**: UCSC / Ensembl에서 선택한 유전자 주변 구간을 바로 받음.
2. **암/정상 서열 패턴 변화 프로토타입**: SEQC-II `SRP162370`에서 체세포 변이 VCF/호출표를 사용하여 동일한 기준 좌표 위에 정상 경로와 변이 반영 경로를 계산함. 처음부터 수십 GB FASTQ를 모두 받는 대신, 먼저 공개 호출표와 관심 유전자 구간으로 분석.
3. **독립 기술 검증**: NIST HG008의 tumor/normal 공개 파일, PacBio Revio / CASTLE long reads.
4. **범암종 탐색**: GDC 공개 체세포 MAF를 사용해 변이 패턴/유전자별 CGR 차이 후보를 찾되 WGS·WES 기원을 구별함.

## 분석 주의

- 정상 **reference DNA**는 개별 환자 정상 조직과 다름.
- SNV 한 개를 적용한 CGR은 종양 전체의 이배체/다배체 서열을 대변하지 않음. allele fraction, purity, CNV, 구조변이 등을 고려해야 함.
- 환자별 암/정상 WGS의 원시 FASTQ에서 실제 서열을 재구성하려면 QC, mapping, calling, haplotype/clone 처리 등 별도의 파이프라인이 필요함.
- 데이터 사용 조건은 각각의 원자료 저장소에서 재확인. 패키지 자체에는 원시 인간 유전체 자료를 포함하지 않음.
