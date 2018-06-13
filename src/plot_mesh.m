% Plot the GD elements in cell array B
function plot_mesh(B)

for k = 1:length(B)
  if isfield(B{k},'isdg')
    continue
  end
  p = B{k}.p;

  x1  = reshape(B{k}.x1, B{k}.Np2, B{k}.Np1);
  x2  = reshape(B{k}.x2, B{k}.Np2, B{k}.Np1);

  if isfield(B{k}, 'E1')
      x1 = B{k}.E2 * x1 * B{k}.E1';
      x2 = B{k}.E2 * x2 * B{k}.E1';
  end
  if isfield(B{k},'isaffine')
    plot(x1(p+1:end-p,[p+1,end-p])' , x2(p+1:end-p,[p+1,end-p])' , 'k');
    hold on
    plot(x1([p+1,end-p],p+1:end-p), x2([p+1,end-p],p+1:end-p), 'k');
    plot(x1([p+1,end-p],[p+1,end-p])' , x2([p+1,end-p],[p+1,end-p])' , 'k', 'LineWidth', 2);
    plot(x1([p+1,end-p],[p+1,end-p]), x2([p+1,end-p],[p+1,end-p]), 'k', 'LineWidth', 2);
  else
    plot(x1(p+1:end-p,p+1:end-p) , x2(p+1:end-p,p+1:end-p) , 'k');
    hold on
    plot(x1(p+1:end-p,p+1:end-p)', x2(p+1:end-p,p+1:end-p)', 'k');
    plot(x1([p+1,end-p],p+1:end-p)' , x2([p+1,end-p],p+1:end-p)' , 'k', 'LineWidth', 2);
    plot(x1(p+1:end-p,[p+1,end-p]), x2(p+1:end-p,[p+1,end-p]), 'k', 'LineWidth', 2);
  end
end

hold off
axis image
