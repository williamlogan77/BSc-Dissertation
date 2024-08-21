import imageio
import os

ele = "Sm"
image_dir = f"./plots/i_process3/{ele}/"

print(len(os.listdir(image_dir)))
images = []
for filename in sorted(os.listdir(image_dir)):
    try:
        # print(f"{image_dir}{filename}")
        images.append(imageio.imread(f"{image_dir}{filename}"))
    except ValueError:
        print(f"BAD {filename}")
imageio.mimsave(f'./plots/i_process3/{ele}/flux.gif', images, duration=0.2)
