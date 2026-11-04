# 공개 WGS 우선 연구자료 실사

> 기존 조사 문서임. 데이터 접근처와 ICGC 포털 종료 등 변경사항은 `OPEN_SEQUENCE_SOURCES_KO.md`를 우선 확인할 것. 공개 *변이 호출표*와 공개 *원시 WGS*는 구분해야 함.

기준일: 2026-11-04. 범위: 특정 암종으로 제한하지 않고, **별도 데이터 접근 승인 없이 사용할 수 있다고 공식 제공처가 확인한 자료만** 선정한다. 최종 사용자별 파일 다운로드 속도·용량·재배포 라이선스는 별도로 확인한다. 통제 접근 원시 WGS를 '공개'로 계산하지 않는다.

## 1. 우선 검토 대상

| 우선 | 자료 | 암종 / 검체 | 공개 WGS 수준 | 메틸화 | 연구 용도 |
|---|---|---|---|---|---|
| 1 | **NIST GIAB HG008** | 췌장암 PDAC 암 세포주, 같은 사람의 정상 췌장·십이지장 조직 | **공개 원시 WGS + 종양/정상 assembly/체세포 변이 벤치마크** | PacBio HiFi 등 direct modification 자료가 공개된 파일에 존재 | 종양 변이 → 서열 CGR 검증 |
| 2 | **CASTLE panel** (PRJNA1086849) | 폐암·유방암 중심 종양/정상 **세포주 6쌍** | **공개 원시 WGS**, Illumina·ONT·PacBio 및 somatic call sets | ONT/PacBio 메틸화 call 공개 | 동일 검체 계열의 서열+메틸화 결합 시험 |
| 3 | **ICGC 25K / PCAWG open tier** | **다암종** ICGC 공개 프로젝트 | **WGS에서 도출된 공개 체세포 변이, CNV/SV 등 결과**; 전체 raw WGS가 아님 | 암종·파일별로 다름, paired methylation 보장 안 됨 | 범암종 서열 패턴 비교 |
| 4 | **SEQC2 HCC1395 / HCC1395BL** (SRP162370) | 유방암 세포주·동일 공여자의 정상 B세포 유래 세포주 | **공개 원시 WGS와 benchmark VCF** | 이 연구의 주목적은 variant benchmarking, 메틸화 통합 코호트 아님 | 체세포 변이 반영 알고리즘 검증 |
| 5 | **COLO829 / COLO829BL** | 흑색종 세포주/정상 B세포 유래 세포주 | **PacBio 공개 HiFi WGS** | 공개 HiFi modification call이 포함된 리드 | 다른 암종 서열 검증 |

CASTLE의 HCC1395는 SEQC2의 HCC1395와 같은 표준 세포주임. **독립 환자 두 명이나 재현 코호트 둘로 세지 않는다.** GEO/SRA의 접근 가능한 파일과 논문 내 생물학적 환자 n 역시 별개이다.

## 2. 확인된 공개 경로

- **NIST GIAB HG008:** https://www.nist.gov/programs-projects/cancer-genome-bottle — NIST 공개 FTP 매니페스트 및 종양/정상 샘플 목록. 기술 설명 논문: https://www.nature.com/articles/s41597-025-05438-2. HG009은 확대 중인 별도 PDAC 전이 세포주 자원으로, 파일 목록 검증 후 우선순위를 정할 것.
- **CASTLE:** https://github.com/CASTLE-Panel/castle ; https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1086849 — 6개 종양/정상 세포주 쌍의 원시 시퀀싱, methylation call 접근 경로.
- **ICGC/PCAWG:** https://docs.icgc-argo.org/docs/data-access/icgc-25k-data — 익명(no-sign-request) S3 공개 객체 저장소. 공개 합의 SNV/indel 결과 예시 경로는 `PCAWG/consensus_snv_indel/`. **EGA/TCGA의 raw WGS는 승인 없이는 안 됨.**
- **SEQC2:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/ ; SRA `SRP162370` — WGS FASTQ, VCF 및 benchmark 공개.
- **COLO829:** https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/COLO829/ — 암/정상 long-read WGS, HiFi modification calls 안내.

### PCAWG 공개 자료 받는 명령 예시

```bash
aws s3 ls s3://icgc25k-open/PCAWG/consensus_snv_indel/ \
  --endpoint-url https://object.genomeinformatics.org --no-sign-request
```

이 경로는 WGS *유래* 변이 결과를 제공한다. 개인별 WGS FASTQ/BAM/CRAM 원시 리드가 공개라는 뜻이 아니다.

## 3. 기준을 충족하지 못해 지금 제외한 주요 임상 코호트

- EAGLE / Sherlock-Lung: WGS와 EPIC 연결 연구이지만 **dbGaP 접근 승인 필요**.
- UCSF WCDT(mCRPC) / DKFZ GBM: WGS+WGBS 등록은 있지만 **EGA 통제 접근**.
- Hartwig metastatic pan-cancer WGS: 임상 원시 및 처리 WGS 모두 공식 자료 접근 절차가 요구됨.
- TCGA-LUAD를 포함한 TCGA raw WGS: **GDC 통제 접근**. 공개 처리된 체세포 변이 파일만 별도 평가 가능.

## 4. 데이터를 앱 연구로 연결하는 원칙

**우선 1: 공개 종양/정상 서열에서 CGR/FCGR 재현성 테스트.** NIST HG008 및 CASTLE의 참조·체세포 변이를 사용하고, SNV-only/indel/SV/CNV 수준을 명시한다. VAF와 copy-number를 무시해 하나의 완벽한 암세포 서열을 만들었다고 주장하지 않는다.

**우선 2: 공개 long-read 메틸화 call을 위치에 맞춰 병합.** CASTLE/NIST 플랫폼별 direct modified-base call의 chemistry 및 motif filtering을 확인한다. 종양/정상 메틸화와 서열이 동일 library/aliquot인지 매니페스트로 검사한다.

**우선 3: 다암종 확장.** PCAWG 오픈 변이 결과를 활용하되 오픈 티어인지, 암종과 검체 정보가 공개되는지, 실제 유전자 좌표가 맞는지 확인한 뒤 병합한다. 관측된 공간 무늬가 염기 빈도와 coverage만 반영하는지 통제 실험 필수.

## 5. 현재 미확정 사항

- 본 세션에서는 각 다운로드 대상의 파일 바이트를 실제 내려받거나 WGS QC를 실행하지 않았음. 공개 카탈로그와 공식 공개 정책을 검증한 단계임.
- 공개 WGS가 있다고 하여 동일 검체의 WGBS도 반드시 존재하는 것은 아님. ONT/PacBio modified base data는 WGBS와 다르며, 플랫폼과 methylation calling 품질이 분석에 영향을 줌.
- ICGC/PCAWG 공개 결과 중 어느 암종에 몇 개 환자의 암·정상 유전체가 모두 사용 가능한지는 프로젝트별 manifest가 필요함. 공식 '전체 암종' 숫자를 임의로 채우지 않음.

이 문서는 **전 세계 모든 코호트의 완전한 목록**이 아니라, 이번 조사에서 공식적으로 공개 접근이 확인된 고우선순위 자료 목록임.
