import zipfile, json, os, time, numpy as np, pandas as pd
from pathlib import Path
base=Path('/mnt/data')
files={'zebrafish':'f2f28d77-0c5c-40aa-aa62-509e40294760.zip','mouse':'5ffeadd4-da71-4a27-a895-13caea41400a.zip'}
report={}
for species,name in files.items():
 t=time.perf_counter();print('Processing',species,flush=True)
 with zipfile.ZipFile(base/name) as z:
  expr_file=next(n for n in z.namelist() if n.endswith('.tsv') and '/metadata_' not in n)
  df=pd.read_csv(z.open(expr_file),sep='\t',index_col=0).astype(np.float32)
  arr=df.to_numpy(dtype=np.float32,copy=False)
  nan=np.isnan(arr).sum()
  finite=np.isfinite(arr).sum()
  negative=(arr<0).sum()
  zero=(arr==0).sum()
  rowmeans=np.nanmean(arr,axis=1)
  rowvar=np.nanvar(np.log1p(np.clip(arr,0,None)),axis=1)
  colmeans=np.nanmean(arr,axis=0)
  colsums=np.nansum(arr,axis=0)
  topids=df.index[np.argsort(rowvar)[-10:][::-1]].tolist()
  report[species]={
   'expression_file':expr_file,
   'samples':int(df.shape[1]),'genes':int(df.shape[0]),
   'unique_gene_ids':int(df.index.nunique()),'non_finite_values':int(arr.size-finite),
   'nan_values':int(nan),'negative_values':int(negative),
   'zero_values':int(zero),'zero_fraction':round(float(zero/arr.size),5),
   'entry_min':round(float(np.nanmin(arr)),6),'entry_max':round(float(np.nanmax(arr)),6),
   'gene_mean_quantiles':dict(zip(['p0','p25','p50','p75','p100'],np.quantile(rowmeans,[0,.25,.5,.75,1]).round(4).tolist())),
   'sample_mean_quantiles':dict(zip(['p0','p25','p50','p75','p100'],np.quantile(colmeans,[0,.25,.5,.75,1]).round(4).tolist())),
   'sample_sum_quantiles':dict(zip(['p0','p25','p50','p75','p100'],np.quantile(colsums,[0,.25,.5,.75,1]).round(2).tolist())),
   'top_10_log1p_variable_gene_ids':topids
  }
  # useful portable top gene list
  variances=pd.DataFrame({'ensembl_gene_id':df.index.to_numpy(),'mean_processed_expr':rowmeans,'log1p_variance':rowvar})
  variances.sort_values('log1p_variance',ascending=False).to_csv(base/f'{species}_gene_variability.csv',index=False)
  pd.DataFrame({'sample_id':df.columns,'mean_processed_expr':colmeans,'sum_processed_expr':colsums}).to_csv(base/f'{species}_sample_summary.csv',index=False)
  del df,arr
 print('Completed',species,round(time.perf_counter()-t,1),'s',report[species],flush=True)
with open(base/'embryo_dataset_qc.json','w') as f:json.dump(report,f,indent=2)
