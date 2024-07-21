import os
import numpy as np
from DPOD.helper import save_obj
from DPOD.pose_block import initial_pose_estimation
from DPOD.create_renderings import create_refinement_inputs
from DPOD.pose_refinement import train_pose_refinement
from DPOD.correspondence_block import train_correspondence_block
from DPOD.create_ground_truth import create_GT_masks, create_UV_XYZ_dictionary, dataset_dir_structure

# Update these paths to your T-LESS dataset and backgrounds
root_dir = "C:\\Users\\ezeun\\sgi2024\\sgi2024\\graph-based-optimal-transport\\datasets\\t-less"
background_dir = "C:\\Users\\ezeun\\sgi2024\\sgi2024\\graph-based-optimal-transport\\datasets\\COCO\\val2017"
intrinsic_matrix = np.array([[572.41140, 0, 325.26110], [0, 573.57043, 242.04899], [0, 0, 1]])

# Update the classes dictionary for T-LESS
classes = {f'{i:02d}': i for i in range(1, 31)}

sensors_train = ['train_primesense', 'train_kinect', 'train_canon']
sensors_test = ['test_primesense', 'test_kinect', 'test_canon']

# Collect all training images
list_all_images = []
for sensor in sensors_train:
    sensor_dir = os.path.join(root_dir, sensor)
    for obj_dir in os.listdir(sensor_dir):
        rgb_path = os.path.join(sensor_dir, obj_dir, 'rgb')
        if os.path.isdir(rgb_path):
            for root, dirs, files in os.walk(rgb_path):
                for file in files:
                    if file.endswith(".jpg"):
                        list_all_images.append(os.path.join(root, file))

num_images = len(list_all_images)
indices = list(range(num_images))
np.random.seed(69)
np.random.shuffle(indices)
split = int(np.floor(0.15 * num_images))
train_idx, test_idx = indices[:split], indices[split:]

save_obj(list_all_images, root_dir + "all_images_adr")
save_obj(train_idx, root_dir + "train_images_indices")
save_obj(test_idx, root_dir + "test_images_indices")

dataset_dir_structure(root_dir)

create_GT_masks(root_dir, background_dir, intrinsic_matrix, classes)
create_UV_XYZ_dictionary(root_dir)

train_correspondence_block(root_dir, classes, epochs=20)

initial_pose_estimation(root_dir, classes, intrinsic_matrix)

create_refinement_inputs(root_dir, classes, intrinsic_matrix)

train_pose_refinement(root_dir, classes, epochs=10)

# Optionally, collect test images for evaluation (not included in training)
list_all_test_images = []
for sensor in sensors_test:
    sensor_dir = os.path.join(root_dir, sensor)
    for scene_dir in os.listdir(sensor_dir):
        rgb_path = os.path.join(sensor_dir, scene_dir, 'rgb')
        if os.path.isdir(rgb_path):
            for root, dirs, files in os.walk(rgb_path):
                for file in files:
                    if file.endswith(".jpg"):
                        list_all_test_images.append(os.path.join(root, file))

save_obj(list_all_test_images, root_dir + "all_test_images_adr")