# %%
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import torch

torch.cuda.empty_cache()
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from utils import AverageMeter
from GIGN import GIGN
from dataset_GIGN import GraphDataset, PLIDataLoader
from config.config_dict import Config
from log.train_logger import TrainLogger
import numpy as np
from utils import *
from sklearn.metrics import mean_squared_error

f = open("loss_MSE_Loss_molecule_lessthan_16A_PDB_Bind_Data.txt", "w")


# %%
def val(model, dataloader, device):
    model.eval()

    pred_list = []
    label_list = []
    for data in dataloader:
        data = data.to(device)
        with torch.no_grad():
            pred = model(data)
            label = data.y

            pred_list.append(pred.detach().cpu().numpy())
            label_list.append(label.detach().cpu().numpy())

    pred = np.concatenate(pred_list, axis=0)
    label = np.concatenate(label_list, axis=0)

    coff = np.corrcoef(pred, label)[0, 1]
    epoch_rmse = np.sqrt(mean_squared_error(label, pred))

    model.train()

    return epoch_rmse, coff


# %%
if __name__ == "__main__":
    cfg = "TrainConfig_GIGN"
    config = Config(cfg)
    args = config.get_config()
    graph_type = args.get("graph_type")
    save_model = args.get("save_model")
    batch_size = 10
    data_root = "/mnt/cephfs/projects/2023110101_Ligand_fitting_to_EM_maps/PDBbind"
    epochs = 2
    repeats = args.get("repeat")
    early_stop_epoch = args.get("early_stop_epoch")

    for repeat in range(repeats):
        args["repeat"] = repeat

        toy_dir = os.path.join(data_root, "PDBBind_Zenodo_6408497")
        toy_df = pd.read_csv(
            os.path.join(toy_dir, "PDB_IDs_with_rdkit_length_less_than_16A.csv")
        ).sample(frac=1.0, random_state=123)
        split_idx = int(0.9 * len(toy_df))
        train_df = toy_df.iloc[:split_idx]
        valid_df = toy_df.iloc[split_idx:]

        train_set = GraphDataset(toy_dir, train_df, graph_type=graph_type, create=False)
        # valid_set = GraphDataset(toy_dir, valid_df, graph_type=graph_type, create=False)

        train_loader = PLIDataLoader(
            train_set, batch_size=batch_size, shuffle=True, drop_last=True
        )
        # valid_loader = PLIDataLoader(valid_set, batch_size=batch_size, shuffle=False, num_workers=4)

        logger = TrainLogger(args, cfg, create=True)
        logger.info(__file__)
        logger.info(f"train data: {len(train_set)}")
        # logger.info(f"valid data: {len(valid_set)}")

        device = torch.device("cuda:0")
        model = GIGN(35, 256).to(device)

        optimizer = optim.Adam(model.parameters(), lr=5e-4, weight_decay=1e-6)
        criterion = nn.MSELoss()
        # epochs = 2

        running_loss = AverageMeter()
        running_acc = AverageMeter()
        running_best_mse = BestMeter("min")
        best_model_list = []

        model.train()
        for epoch in range(epochs):
            for data in train_loader:
                data = data.to(device)
                print(data.size())
                pred = model(data)
                label = data.y
                print("PRED", pred.size())
                print("LABEL", label.size())

                loss = criterion(pred, label)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                running_loss.update(loss.item(), label.size(0))

            epoch_loss = running_loss.get_average()
            epoch_rmse = np.sqrt(epoch_loss)
            running_loss.reset()
            print(epoch_loss)
            f.write(str(epoch_loss) + "\n")

            # start validating
            # valid_rmse, valid_pr = val(model, valid_loader, device)
            # msg = "epoch-%d, train_loss-%.4f, train_rmse-%.4f, valid_rmse-%.4f, valid_pr-%.4f" \
            #        % (epoch, epoch_loss, epoch_rmse, valid_rmse, valid_pr)
            # logger.info(msg)

            # if valid_rmse < running_best_mse.get_best():
            #    running_best_mse.update(valid_rmse)
            #    if save_model:
            #        msg = "epoch-%d, train_loss-%.4f, train_rmse-%.4f, valid_rmse-%.4f, valid_pr-%.4f" \
            #        % (epoch, epoch_loss, epoch_rmse, valid_rmse, valid_pr)
            #        model_path = os.path.join(logger.get_model_dir(), msg + '.pt')
            #        best_model_list.append(model_path)
            #        save_model_dict(model, logger.get_model_dir(), msg)
            # else:
            #    count = running_best_mse.counter()
            #    if count > early_stop_epoch:
            #        best_mse = running_best_mse.get_best()
            #        msg = "best_rmse: %.4f" % best_mse
            #        logger.info(f"early stop in epoch {epoch}")
            #        logger.info(msg)
            #        break_flag = True
            #        break

f.close()
import joblib

joblib.dump(model, "model_saved_molecule_lessthan_16A_PDB_Bind_Data.pkl")

# %%
