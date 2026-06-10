# Colab pathway pipeline

Run these cells in order in Google Colab.

```python
%cd /content
!rm -rf istar-pathway-loss
!git clone https://github.com/verysleepyguy/istar-pathway-loss.git
%cd istar-pathway-loss
```

```python
!pip install -r requirements.txt
!pip install scanpy squidpy anndata
```

```python
import torch
print("CUDA available:", torch.cuda.is_available())
```

```python
!python make_squidpy_demo_data.py --prefix data/visium_mouse_brain/
!python validate_inputs.py data/visium_mouse_brain/ --pathways data/pathways-demo.tsv
!ls -lh data/visium_mouse_brain/
```

```python
!bash download_checkpoints.sh
```

Use a short smoke test first:

```python
!EPOCHS=5 N_STATES=1 DEVICE=cuda bash run_pathway.sh data/visium_mouse_brain/ data/pathways-demo.tsv
```

If the smoke test finishes, run the full pathway model:

```python
!rm -rf data/visium_mouse_brain/states-pathway data/visium_mouse_brain/cnts-pathway-super
!EPOCHS=400 N_STATES=5 DEVICE=cuda bash run_pathway.sh data/visium_mouse_brain/ data/pathways-demo.tsv
```

Check outputs:

```python
!ls -lh data/visium_mouse_brain/cnts-pathway-super/
!ls -lh data/visium_mouse_brain/cnts-pathway-super-plots/
!cat data/visium_mouse_brain/pathway-names.txt
!head data/visium_mouse_brain/pathway-matched-genes.tsv
```
