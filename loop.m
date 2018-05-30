
w=100;                              %检测区域的宽                              %静态节点数
l=100;
n=length(x_x)

[v,c]=voronoin([x_x;y_y]');
%v就是所有的voronoi边的顶点坐标，c{i}就是第i个点周围的voronoi边的顶点下标

axis([0,l,0,w])
axis equal
%比如说第1个点(0,0)周围的voronoi边的顶点就是横坐标：v(c{1},1)，纵坐标：v(c{1},2)
%对于坐标为(inf,inf)的点就是无穷远点.
k_slove=[];
b_slove=[];
x_node_temp=[];
y_node_temp=[];
for n=1:n
a=[v(c{n},1),v(c{n},2)];
    for i=1:1:length(v(c{n}))               %循环每个静态节点形成的泰森多边形的各个顶点求出坐标
        if(i<length(v(c{n})))
            x_x0=a(i,1);
            y_y0=a(i,2);
            x_x1=a(i+1,1);
            y_y1=a(i+1,2);
                p1 = [x_x0 y_y0];
                p2 = [x_x1 y_y1];
                p3=[x_x(n) y_y(n)];
                l3=norm(p1-p2);                 %泰森多边形一个边的长度
                l1=norm(p3-p1);                 %静态节点到泰森多边形边的其中一端的距离
                l2=norm(p3-p2);                 %静态节点到泰森多边形边的另一端的距离
                k=(y_y1-y_y0)/(x_x1-x_x0);              %求出直线方程的K
                b=y_y0-k*x_x0;                                            %求出直线方程的b
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %此处应判断k_slove和b_slove中是否已经有对应的k和b，如果没有则进行下一步，如果有则跳过本次循环
            flog=0;
            for p=1:1:length(k_slove)                
                if(k_slove(p)==k&&b_slove(p)==b)
                    flog=1;
                end 
                if(l3<1.5)
                    flog=1;
                end 
            end
            if(flog==1)
                disp('这句话被执行了')
                continue
                
            end

            k_slove=[k_slove k];                                    %存储循环过程中产生的k
            b_slove=[b_slove b];                                    %存储循环过程中产生的b
           
            d=abs((k*x_x(n)-y_y(n)+b)/sqrt(k*k+1));                     %点到直线的距离
            if(r<d)                                                 %判断两邻居节点的距离是否大于二倍的Rs
                al=acos(abs((l1*l1+l2*l2-l3*l3))/2*l1*l2);          %漏洞探测  计算圆心角
                fugailv=al*r*r/l3*d;                                %漏洞探测  计算覆盖率
                 disp(fugailv)
                    if(fugailv<0.9)                                 %漏洞探测  覆盖率小于百分之八十即修补漏洞
                        a0=x_x(n);
                        a1=y_y(n);
                        b0=x_x0;
                        b1=y_y0;
                        c0=x_x1;
                        c1=y_y1;
                        Dab=sqrt((b0-a0)^2+(b1-a1)^2);      %计算角平分线与对边的交点
                        Dac=sqrt((c0-a0)^2+(c1-a1)^2);      %计算角平分线与对边的交点
                        value=Dab/Dac;                                  %计算角平分线与对边的交点
                        x_node=(b0+value*c0)/(1+value);                 %计算角平分线与对边的交点的横坐标
                        y_node=(b1+value*c1)/(1+value);                 %计算角平分线与对边的交点的纵坐标
                        
                        fleg=0;
                        for q=1:1:length(x_node_temp)                
                            if(x_node_temp(q)==x_node||y_node_temp(q)==y_node)
                                fleg=1;
                            end 
                        end
                        if(fleg==1)
                            continue
                        end
                        
                         x_node_temp=[x_node_temp x_node];                                  
                        y_node_temp=[y_node_temp y_node];
                       
                            fill(x_node+r*cos(theta),y_node+r*sin(theta),'r','EdgeColor','w','facealpha',0.2);%画圆
                            hold on
                            plot(x_node,y_node,'k.');           %画静态节点
                            hold on

                        disp('这句话没有执行到')
                    end                                             %漏洞探测
            else
                    disp('此处R>d')   
            end                                                     %漏洞探测
        end
    end
end
axis([0,l,0,w])
for  i=1:1:n                        %for循环画出静态节点的感知圆
    fill(x_x(i)+r*cos(theta),y_y(i)+r*sin(theta),'r','EdgeColor','w','facealpha',0.2);%画圆
    hold on
    plot(x_x(i),y_y(i),'b.');
    voronoi(x_x,y_y);axis([0 l 0 w])
hold on
  figure(3)
  hold on
end
axis equal
axis([0,l,0,w])