# 配置 Galaxy

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

```bash
# 配置主设置文件
(galaxy)$ cd galaxy
(galaxy)$ cp galaxy.yml.sample galaxy.yml
```

**建立管理员**

注册一个新用户，将其email地址作为管理员用户添加到`galaxy.yml`中。

```bash
# 添加admin，将实际的注册email信息填入，如果要多个管理员，用逗号区分
admin_users: abc@123.com,myuser@new.com
```

配置文件中添加管理员后，在主页菜单栏中会出现`管理员`菜单，点击进入。
