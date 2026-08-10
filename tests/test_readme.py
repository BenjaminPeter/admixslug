def test_quickstart_bam(script_runner):
    cmd = "admixslug-bam --bam tests/data/FQ.L30MQ25_bigsteffi52.bam --out tests/res/FQ_bigsteffi52.in.xz  --ref tests/data/ref_bigsteffi.csv.xz"
    args = cmd.split()
    ret = script_runner.run(args, cwd="./")
    assert ret.success

def test_quickstart_slug(script_runner):
    cmd = """admixslug --infile tests/data/FQ_bigsteffi52.in.xz 
        --ref tests/data/ref_bigsteffi.csv.xz 
        -o tests/res/FQ_bigsteffi52  
        --states ALT VIN CHA DEN  
        --cont-id EUR  
        --ancestral PAN 
        --ll-tol 0.01  
        --ptol 0.001   
        --max-iter 100 
        --filter-pos 50 
        --filter-ancestral  
        --len-bin-size 2000
        --jk-resamples 10
        --output-jk-sfs
        --output-fstats"""
    args = cmd.split()
    ret = script_runner.run(args, cwd="./")
    assert ret.success
