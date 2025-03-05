###
Create and run the docker file:
```
docker build --no-cache -t ext .

docker run -it ext
```
Run the experiments:

```
source activate quark_install

python instance.py
```

Create Plots:
```
Rscript plots.r
```