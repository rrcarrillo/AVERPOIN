% num_vs_an_plot(C,max_points) plots a graph showing the fitness value of
% C calculated analyticaly and numerically
% C is a matrix specifying the input representation and max_points is the
% maximum number of point (rectangles) to evaluate for the numerical
% integration
function num_vs_an_plot(C,max_points)

int_an=int_res(C);

num_fitness_list=[];
x_coords=1:1:max_points;
for num_points=x_coords
    int_num=int_res_num(C,num_points,false);
    num_fitness_list=[num_fitness_list int_num.norm_integ];
end
plot(x_coords,num_fitness_list,'x-r','LineWidth',2)
hold on
plot(x_coords,ones(1,length(x_coords))*int_an.norm_integ,'b','LineWidth',2)
%title('Numerical vs analytical calculation comparison')
xlabel('Num. resolution N (points/dim)')
ylabel('Normalized representation error')
legend('Numerical integration Ir\_numN(C,N)','Analytical calculation IrN(C)','Location','SouthEast')
axis([x_coords(1) x_coords(end) 0 0.03])
