
w=100;                              %�������Ŀ�                              %��̬�ڵ���
l=100;
n=length(x_x)

[v,c]=voronoin([x_x;y_y]');
%v�������е�voronoi�ߵĶ������꣬c{i}���ǵ�i������Χ��voronoi�ߵĶ����±�

axis([0,l,0,w])
axis equal
%����˵��1����(0,0)��Χ��voronoi�ߵĶ�����Ǻ����꣺v(c{1},1)�������꣺v(c{1},2)
%��������Ϊ(inf,inf)�ĵ��������Զ��.
k_slove=[];
b_slove=[];
x_node_temp=[];
y_node_temp=[];
for n=1:n
a=[v(c{n},1),v(c{n},2)];
    for i=1:1:length(v(c{n}))               %ѭ��ÿ����̬�ڵ��γɵ�̩ɭ����εĸ��������������
        if(i<length(v(c{n})))
            x_x0=a(i,1);
            y_y0=a(i,2);
            x_x1=a(i+1,1);
            y_y1=a(i+1,2);
                p1 = [x_x0 y_y0];
                p2 = [x_x1 y_y1];
                p3=[x_x(n) y_y(n)];
                l3=norm(p1-p2);                 %̩ɭ�����һ���ߵĳ���
                l1=norm(p3-p1);                 %��̬�ڵ㵽̩ɭ����αߵ�����һ�˵ľ���
                l2=norm(p3-p2);                 %��̬�ڵ㵽̩ɭ����αߵ���һ�˵ľ���
                k=(y_y1-y_y0)/(x_x1-x_x0);              %���ֱ�߷��̵�K
                b=y_y0-k*x_x0;                                            %���ֱ�߷��̵�b
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %�˴�Ӧ�ж�k_slove��b_slove���Ƿ��Ѿ��ж�Ӧ��k��b�����û���������һ�������������������ѭ��
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
                disp('��仰��ִ����')
                continue
                
            end

            k_slove=[k_slove k];                                    %�洢ѭ�������в�����k
            b_slove=[b_slove b];                                    %�洢ѭ�������в�����b
           
            d=abs((k*x_x(n)-y_y(n)+b)/sqrt(k*k+1));                     %�㵽ֱ�ߵľ���
            if(r<d)                                                 %�ж����ھӽڵ�ľ����Ƿ���ڶ�����Rs
                al=acos(abs((l1*l1+l2*l2-l3*l3))/2*l1*l2);          %©��̽��  ����Բ�Ľ�
                fugailv=al*r*r/l3*d;                                %©��̽��  ���㸲����
                 disp(fugailv)
                    if(fugailv<0.9)                                 %©��̽��  ������С�ڰٷ�֮��ʮ���޲�©��
                        a0=x_x(n);
                        a1=y_y(n);
                        b0=x_x0;
                        b1=y_y0;
                        c0=x_x1;
                        c1=y_y1;
                        Dab=sqrt((b0-a0)^2+(b1-a1)^2);      %�����ƽ������ԱߵĽ���
                        Dac=sqrt((c0-a0)^2+(c1-a1)^2);      %�����ƽ������ԱߵĽ���
                        value=Dab/Dac;                                  %�����ƽ������ԱߵĽ���
                        x_node=(b0+value*c0)/(1+value);                 %�����ƽ������ԱߵĽ���ĺ�����
                        y_node=(b1+value*c1)/(1+value);                 %�����ƽ������ԱߵĽ����������
                        
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
                       
                            fill(x_node+r*cos(theta),y_node+r*sin(theta),'r','EdgeColor','w','facealpha',0.2);%��Բ
                            hold on
                            plot(x_node,y_node,'k.');           %����̬�ڵ�
                            hold on

                        disp('��仰û��ִ�е�')
                    end                                             %©��̽��
            else
                    disp('�˴�R>d')   
            end                                                     %©��̽��
        end
    end
end
axis([0,l,0,w])
for  i=1:1:n                        %forѭ��������̬�ڵ�ĸ�֪Բ
    fill(x_x(i)+r*cos(theta),y_y(i)+r*sin(theta),'r','EdgeColor','w','facealpha',0.2);%��Բ
    hold on
    plot(x_x(i),y_y(i),'b.');
    voronoi(x_x,y_y);axis([0 l 0 w])
hold on
  figure(3)
  hold on
end
axis equal
axis([0,l,0,w])