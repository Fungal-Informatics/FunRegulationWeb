import { MigrationDefinitionWithName } from "umzug/lib/migrationsList";
import { zap } from "../helpers/sql";

const migration: MigrationDefinitionWithName = {
	name: "05-bookmarks",
	async up() {
		await zap.transaction.single/*sql*/`
            create table 
                bookmarks
            ( 
                id uuid primary key,
                title text not null,
                description text,
                url text not null,
				"boardId" uuid references boards(id) not null,
                "createdAt" timestamptz not null,
                "updatedAt" timestamptz not null,
				"isArchived" boolean not null
            );
        `;
	},
	async down() {
		await zap.transaction.single/*sql*/`
            drop table bookmarks
        `;
	},
};

export default migration;
