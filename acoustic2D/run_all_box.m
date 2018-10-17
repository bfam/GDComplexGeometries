diary('run_all_box.diary')
clear
for p = 1:5
  clearvars -except p
  GDBox_curved
end
clear
for p = 1:5
  clearvars -except p
  GDBox_curved_extrapolation
end
diary off
