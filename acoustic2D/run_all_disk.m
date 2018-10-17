diary('run_all_disk.diary')
clear
for p = 1:5
  clearvars -except p
  DGDisk2D
end
clear
for p = 1:5
  clearvars -except p
  GDDGDisk2D
end
clear
for p = 1:5
  clearvars -except p
  GDDGDisk2D_extrapolation
end
clear
for p = 1:5
  clearvars -except p
  GDDisk2D
end
clear
for p = 1:5
  clearvars -except p
  GDDisk2D_extrapolation
end
diary off
